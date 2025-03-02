import pandas as pd
import numpy as np
from typing import Dict, List, Union, Tuple
import re

class VariantPrioritizer:
    def __init__(self):
        # Scoring weights with higher emphasis on clinical evidence
        self.weights = {
            'functional_impact': 0.15,    
            'population_freq': 0.15,      # Increased weight for rare variants
            'clinical_evidence': 0.40,    # Highest weight for clinical evidence
            'computational': 0.15,        
            'conservation': 0.06,         # Slightly reduced
            'splicing': 0.10             # Slightly reduced
        }
        
        # Define impact categories for functional regions
        self.functional_impact = {
            'HIGH': {
                'location': ['exonic', 'exonic;splicing', 'splicing'],
                'effect': ['frameshift deletion', 'frameshift insertion', 'stopgain', 
                          'stoploss', 'startloss', 'startgain', 'nonsynonymous SNV']
            },
            'MODERATE': {
                'location': ['UTR5', 'UTR3', 'UTR5;UTR3'],
                'effect': ['nonframeshift deletion', 'nonframeshift insertion']
            },
            'LOW': {
                'location': ['intronic', 'ncRNA_exonic', 'ncRNA_splicing'],
                'effect': ['synonymous SNV']
            },
            'MINIMAL': {
                'location': ['upstream', 'downstream', 'intergenic', 'ncRNA_intronic'],
                'effect': ['unknown']
            }
        }

        # Define column names (with exact spacing)
        self.column_names = {
            'clinvar': 'clinvar: Clinvar ',
            'intervar': ' InterVar: InterVar and Evidence ',
            'freq_gnomad': 'Freq_gnomAD_genome_ALL'
        }
        
        # Important gene/disease patterns - expanded list focusing on actionable genes
        self.high_impact_gene_patterns = [
            # Cancer-related terms
            'leukemia', 'cancer', 'tumor', 'carcinom', 'lymphom', 'sarcom', 
            'melanom', 'blastom', 'neoplasia', 'dysplasia', 'oncogene', 'metasta',
            
            # Critical gene functions
            'tumor suppressor', 'dna repair', 'cell cycle', 'checkpoint',
            
            # Specific major cancer genes
            'brca', 'tp53', 'egfr', 'kras', 'nras', 'braf', 'pik3ca', 'apc', 'rb1',
            'bcr', 'abl', 'bcr-abl', 'myc', 'erbb2', 'her2', 'ret', 'alk', 'chek2',
            'palb2', 'atm', 'msh2', 'mlh1', 'pms2', 'epcam', 'kit', 'pdgfra',
            
            # Critical disease terms
            'cardiomyopathy', 'channelopath', 'sudden death', 'arrhythmia', 
            'long qt', 'hypertrophic', 'dilated', 'restrictive',
            
            # Neurodevelopmental/degenerative
            'alzheimer', 'parkinson', 'huntington', 'ataxia', 'spastic', 
            'intellectual disability', 'autism', 'epilepsy'
        ]

    def _is_high_impact_variant(self, row: pd.Series) -> bool:
        """Determine if variant has high molecular impact"""
        high_impact_effects = {
            'frameshift', 'nonsense', 'splice-site', 'start-loss',
            'stop-gain', 'stop-loss', 'startgain'
        }
        effect = str(row['ExonicFunc.refGene']).lower()
        
        # Check ExonicFunc.refGene for frameshift and other high-impact variants
        if any(impact in effect for impact in high_impact_effects):
            return True
            
        # Check Func.refGene for splicing variants
        func = str(row['Func.refGene']).lower()
        if 'splicing' in func:
            return True
            
        return False

    def _has_disease_annotation(self, row: pd.Series) -> bool:
        """Check if variant is in a known disease gene"""
        disease_fields = ['Disease_description', 'OMIM', 'Orphanet', 'Gene_full_name', 'CLNDN', 'Otherinfo']
        for field in disease_fields:
            if field in row and pd.notna(row[field]) and row[field] != '.':
                return True
        return False
        
    def _is_in_high_impact_gene(self, row: pd.Series) -> bool:
        """Check if variant is in a high-impact gene (cancer, cardiac, etc.)"""
        # Fields to check for high-impact patterns
        gene_info_fields = ['Ref.Gene', 'Disease_description', 'OMIM', 'Orphanet', 
                          'Gene_full_name', 'CLNDN', 'Otherinfo', 'Phenotype_MIM']
        
        # Compile all text to search
        text_to_check = ""
        for field in gene_info_fields:
            if field in row and pd.notna(row[field]) and row[field] != '.':
                text_to_check += str(row[field]).lower() + " "
        
        # Check for any high-impact pattern
        for pattern in self.high_impact_gene_patterns:
            if pattern.lower() in text_to_check:
                return True
                
        return False
        
    def _has_pvs1_evidence(self, row: pd.Series) -> bool:
        """Check if variant has PVS1 evidence in InterVar"""
        intervar_str = str(row[self.column_names['intervar']]).lower()
        return "pvs1=1" in intervar_str

    def score_functional_impact(self, row: pd.Series) -> float:
        """Score variant based on functional impact"""
        func = str(row['Func.refGene'])
        effect = str(row['ExonicFunc.refGene'])
        
        # Initialize score
        score = 0.0
        
        # Score based on location
        for impact, details in self.functional_impact.items():
            if any(loc in func for loc in details['location']):
                if impact == 'HIGH':
                    score = 1.0
                elif impact == 'MODERATE':
                    score = 0.7
                elif impact == 'LOW':
                    score = 0.4
                else:  # MINIMAL
                    score = 0.1
                break
        
        # Adjust score based on effect if present
        if pd.notna(effect) and effect != '.':
            effect_score = 0.0
            for impact, details in self.functional_impact.items():
                if effect in details['effect']:
                    if impact == 'HIGH':
                        effect_score = 1.0
                    elif impact == 'MODERATE':
                        effect_score = 0.7
                    elif impact == 'LOW':
                        effect_score = 0.4
                    break
            # Combine location and effect scores
            score = max(score, effect_score)
        
        return score

    def score_population_frequency(self, freq: Union[float, str], row: pd.Series = None) -> float:
        """Context-aware population frequency scoring with gene and variant importance considerations"""
        if pd.isna(freq) or freq == '.':
            return 0.7  # Prioritize investigation of unknown variants
        
        try:
            freq = float(freq)
            
            # Check for special cases that need frequency penalty adjustment
            has_pvs1 = False
            in_high_impact_gene = False
            is_high_impact_variant = False
            
            if row is not None:
                has_pvs1 = self._has_pvs1_evidence(row)
                in_high_impact_gene = self._is_in_high_impact_gene(row)
                is_high_impact_variant = self._is_high_impact_variant(row)
            
            # Standard scoring scale for most variants
            if freq == 0:
                return 1.0  # Novel variant
            elif freq <= 0.0001:  # Ultra-rare: 0.01%
                return 0.95
            elif freq <= 0.0005:  # Very rare: 0.05%
                return 0.90
            elif freq <= 0.001:   # Rare: 0.1%
                return 0.85
            elif freq <= 0.005:   # Uncommon: 0.5%
                return 0.70
                
            # Critical PVS1 variant in high-impact gene gets reduced frequency penalty
            if has_pvs1 and in_high_impact_gene and is_high_impact_variant:
                if freq <= 0.01:    # 1%
                    return 0.65     # Much reduced penalty
                elif freq <= 0.05:  # 5%
                    return 0.55     # Significantly reduced penalty 
                else:
                    return 0.45     # Still prioritize over truly common variants
            
            # PVS1 variants get somewhat reduced frequency penalty
            elif has_pvs1 and is_high_impact_variant:
                if freq <= 0.01:    # 1%
                    return 0.55     # Reduced penalty
                elif freq <= 0.05:  # 5%
                    return 0.40     # Less reduced penalty
                else:
                    return 0.30     # Still give some priority
            
            # Regular variants follow standard frequency penalties
            elif freq <= 0.01:    # 1%
                return 0.40
            elif freq <= 0.05:    # 5%
                return 0.20
            else:
                return 0.10
                
        except ValueError:
            return 0.7

    def parse_clinvar(self, clinvar_str: str) -> float:
        """Comprehensive ClinVar parsing for all possible classification scenarios"""
        if pd.isna(clinvar_str) or clinvar_str == '.' or clinvar_str == 'CLNSIG' or clinvar_str == 'UNK':
            return 0.5

        # Remove prefixes and standardize text
        clinvar_str = clinvar_str.lower()
        for prefix in ['clinvar: ', 'clinvar ']:
            if clinvar_str.startswith(prefix):
                clinvar_str = clinvar_str[len(prefix):]
                break
                
        # Handle encoding variations
        clinvar_str = clinvar_str.replace('\\x2c', ',')
        clinvar_str = clinvar_str.replace(',_', '_')
        
        # Primary classification scores
        primary_scores = {
            # Pathogenic variants (highest priority)
            'pathogenic': 1.0,
            'pathogenic/likely_pathogenic': 0.95,
            'pathogenic/likely_pathogenic/likely_risk_allele': 0.93,
            'pathogenic/likely_risk_allele': 0.90,
            'pathogenic/pathogenic_low_penetrance': 0.90,
            'pathogenic/likely_pathogenic/pathogenic_low_penetrance': 0.90,
            'pathogenic/likely_pathogenic/pathogenic_low_penetrance/established_risk_allele': 0.92,
            'pathogenic_low_penetrance': 0.85,
            
            # Likely pathogenic variants
            'likely_pathogenic': 0.85,
            'likely_pathogenic/likely_risk_allele': 0.83,
            'likely_pathogenic_low_penetrance': 0.80,
            'likely_risk_allele': 0.70,
            
            # Uncertain and conflicting variants
            'uncertain_significance': 0.50,
            'uncertain_significance/uncertain_risk_allele': 0.50,
            'uncertain_risk_allele': 0.50,
            'conflicting_classifications_of_pathogenicity': 0.50,
            'association_not_found': 0.50,
            'no_classifications_from_unflagged_records': 0.50,
            'not_provided': 0.50,
            
            # Likely benign variants
            'likely_benign': 0.20,
            
            # Benign variants (lowest priority)
            'benign/likely_benign': 0.15,
            'benign': 0.10
        }
        
        # Modifier impacts (adjust the base score)
        modifier_impact = {
            'risk_factor': 0.05,           # Slightly increase for risk factors
            'affects': 0.03,               # Minor increase for variants that affect phenotype
            'association': 0.02,           # Small increase for associated variants
            'established_risk_allele': 0.04, # Notable increase for established risk
            'confers_sensitivity': 0.01,   # Minimal increase for sensitivity
            'other': 0.00,                 # No change for "other" annotation
            'drug_response': -0.02,        # Minor decrease for drug response
            'protective': -0.05            # More significant decrease for protective variants
        }
        
        # Split into main classification and modifiers
        parts = clinvar_str.split('|')
        main_classification = parts[0]
        modifiers = parts[1:] if len(parts) > 1 else []
        
        # Initialize base score
        base_score = 0.5  # Default to VUS if no match found
        
        # Try for exact match with primary classification
        if main_classification in primary_scores:
            base_score = primary_scores[main_classification]
        else:
            # Try for partial matches, prioritizing longer matches first
            sorted_classifications = sorted(primary_scores.keys(), key=len, reverse=True)
            matched = False
            
            for classification in sorted_classifications:
                if classification in main_classification:
                    base_score = primary_scores[classification]
                    matched = True
                    break
            
            # Special handling for combined classifications not explicitly listed
            if not matched:
                if 'pathogenic' in main_classification:
                    if 'likely_pathogenic' in main_classification:
                        base_score = 0.95  # Combined pathogenic/likely_pathogenic
                    elif 'low_penetrance' in main_classification:
                        base_score = 0.85  # Pathogenic with low penetrance
                    else:
                        base_score = 0.90  # Generic pathogenic combination
                elif 'likely_pathogenic' in main_classification:
                    base_score = 0.85  # Generic likely pathogenic combination
                elif 'benign' in main_classification and 'likely_benign' in main_classification:
                    base_score = 0.15  # Combined benign/likely_benign
                elif 'benign' in main_classification:
                    base_score = 0.10  # Generic benign combination
                elif 'likely_benign' in main_classification:
                    base_score = 0.20  # Generic likely benign combination
                
        # Calculate modifier adjustment
        modifier_adjustment = 0.0
        
        # Check main classification for embedded modifiers
        for modifier, adjustment in modifier_impact.items():
            if modifier in main_classification and 'pathogenic' not in modifier and 'benign' not in modifier:
                modifier_adjustment += adjustment
        
        # Check explicit modifiers (parts after the pipe)
        for modifier_part in modifiers:
            for modifier, adjustment in modifier_impact.items():
                if modifier in modifier_part:
                    modifier_adjustment += adjustment
        
        # Stand-alone classifications that appear without primary classification
        standalone_modifiers = {
            'affects': 0.60,
            'association': 0.55,
            'drug_response': 0.45,
            'risk_factor': 0.65,
            'protective': 0.40,
            'confers_sensitivity': 0.50,
            'other': 0.50
        }
        
        # If the main part is just a modifier with no primary classification
        if base_score == 0.5 and main_classification in standalone_modifiers:
            base_score = standalone_modifiers[main_classification]
        
        # Ensure final score is between 0 and 1
        final_score = max(0.0, min(1.0, base_score + modifier_adjustment))
        
        return final_score

    def _count_evidence(self, intervar_str: str, evidence_type: str) -> int:
        """Count specific evidence type occurrences"""
        pattern = f"{evidence_type}=\\[([^\\]]+)\\]"
        match = re.search(pattern, intervar_str)
        if not match:
            return 0
        return sum(1 for x in match.group(1).split(',') if x.strip() == '1')

    def parse_intervar(self, intervar_str: str) -> float:
        """Enhanced InterVar parsing with better evidence integration"""
        if pd.isna(intervar_str) or intervar_str == '.':
            return 0.5

        intervar_str = intervar_str.lower()
        
        # Base classification scores
        classification_scores = {
            'pathogenic': 0.95,
            'likely pathogenic': 0.85,
            'uncertain significance': 0.5,
            'likely benign': 0.2,
            'benign': 0.1
        }
        
        # Initialize scores
        base_score = 0.5
        evidence_score = 0.0
        
        # Get base classification score
        for classification, score in classification_scores.items():
            if classification in intervar_str:
                base_score = score
                break
        
        # Parse evidence
        evidence = {
            'pvs': self._count_evidence(intervar_str, 'pvs'),
            'ps': self._count_evidence(intervar_str, 'ps'),
            'pm': self._count_evidence(intervar_str, 'pm'),
            'pp': self._count_evidence(intervar_str, 'pp'),
            'ba': self._count_evidence(intervar_str, 'ba'),
            'bs': self._count_evidence(intervar_str, 'bs'),
            'bp': self._count_evidence(intervar_str, 'bp')
        }
        
        # Calculate evidence adjustments - ENHANCED to give PVS1 standalone strong weight
        if evidence['pvs'] > 0:
            # Start with a significant boost for PVS1 alone
            evidence_score = 0.15
            
            # Add additional boost with other evidence
            if evidence['ps'] >= 1:
                evidence_score += 0.05  # More with strong evidence
            
            if evidence['pm'] >= 2:
                evidence_score += 0.05  # More with multiple moderate evidence
            elif evidence['pm'] >= 1:
                evidence_score += 0.03  # Slight boost with one moderate
                
            if evidence['pp'] >= 2:
                evidence_score += 0.02  # Minor boost with supporting evidence
        
        # Other evidence combinations (without PVS1)
        elif evidence['ps'] >= 2:
            evidence_score = 0.15
        elif evidence['ps'] == 1 and evidence['pm'] >= 3:
            evidence_score = 0.12
        elif evidence['ps'] == 1 and evidence['pm'] >= 2:
            evidence_score = 0.10
        elif evidence['ps'] == 1:
            evidence_score = 0.08
        
        # Strong benign evidence
        if evidence['ba'] > 0 or evidence['bs'] >= 2:
            evidence_score -= 0.20
        elif evidence['bs'] == 1 and evidence['bp'] >= 1:
            evidence_score -= 0.15
        elif evidence['bs'] == 1 and evidence['pvs'] == 0 and evidence['ps'] == 0:
            # If there's one benign criterion and no strong pathogenic evidence,
            # reduce score more significantly
            evidence_score -= 0.12
        
        return min(1.0, max(0.0, base_score + evidence_score))

    def score_clinical_evidence(self, row: pd.Series) -> float:
        """Enhanced clinical evidence scoring with better handling of PVS1"""
        clinvar_score = self.parse_clinvar(row[self.column_names['clinvar']])
        intervar_score = self.parse_intervar(row[self.column_names['intervar']])
        
        # Check variant characteristics
        is_high_impact = self._is_high_impact_variant(row)
        has_disease_annotation = self._has_disease_annotation(row)
        has_pvs1 = self._has_pvs1_evidence(row)
        in_high_impact_gene = self._is_in_high_impact_gene(row)
        
        # Special handling for PVS1 variants
        if has_pvs1:
            # For critical genes with PVS1, heavily prioritize InterVar score
            if in_high_impact_gene:
                intervar_weight = 0.85
                clinvar_weight = 0.15
                
                # For uncertain/conflicting ClinVar with PVS1 in important genes, trust InterVar even more
                if 'uncertain' in str(row[self.column_names['clinvar']]).lower() or 'conflict' in str(row[self.column_names['clinvar']]).lower():
                    intervar_weight = 0.90
                    clinvar_weight = 0.10
                
                return (intervar_score * intervar_weight) + (clinvar_score * clinvar_weight)
            
            # For other PVS1 variants, still prioritize InterVar but less extremely
            else:
                intervar_weight = 0.75
                clinvar_weight = 0.25
                return (intervar_score * intervar_weight) + (clinvar_score * clinvar_weight)
        
        # Adjust weighting based on evidence strength for non-PVS1 variants
        if is_high_impact and has_disease_annotation:
            if clinvar_score >= 0.85:  # Strong ClinVar evidence
                return (intervar_score * 0.3) + (clinvar_score * 0.7)  # More weight to ClinVar
            elif intervar_score >= 0.85:  # Strong InterVar evidence
                return (intervar_score * 0.7) + (clinvar_score * 0.3)  # More weight to InterVar
        
        # Default weighting
        return (intervar_score * 0.6) + (clinvar_score * 0.4)  # Slightly increased InterVar influence

    def score_computational_predictions(self, row: pd.Series) -> float:
        """Enhanced computational prediction scoring"""
        scores = []
        weights = []
        
        # CADD score (higher confidence for extreme values)
        if pd.notna(row['CADD_phred']) and row['CADD_phred'] != '.':
            cadd = float(row['CADD_phred'])
            score = min(1.0, cadd / 40)
            confidence = min(1.0, cadd / 50)
            scores.append(score)
            weights.append(1.2 * confidence)
        
        # REVEL score
        if pd.notna(row['REVEL_score']) and row['REVEL_score'] != '.':
            revel = float(row['REVEL_score'])
            scores.append(revel)
            weights.append(1.1)
        
        # Other predictors
        standard_predictors = {
            'SIFT_score': lambda x: 1 - float(x),
            'AlphaMissense_pred': lambda x: 1.0 if 'pathogenic' in str(x).lower() else 0.2,
            'ClinPred_score': lambda x: float(x)
        }
        
        for pred, score_func in standard_predictors.items():
            if pd.notna(row[pred]) and row[pred] != '.':
                try:
                    scores.append(score_func(row[pred]))
                    weights.append(1.0)
                except (ValueError, TypeError):
                    continue
        
        return np.average(scores, weights=weights) if scores else 0.5

    def score_conservation(self, row: pd.Series) -> float:
        """Score variant based on conservation metrics"""
        scores = []
        
        if pd.notna(row['phyloP100way_vertebrate']) and row['phyloP100way_vertebrate'] != '.':
            phylop = float(row['phyloP100way_vertebrate'])
            scores.append(min(1.0, max(0, (phylop + 5) / 10)))
        
        if pd.notna(row['GERP++_RS']) and row['GERP++_RS'] != '.':
            gerp = float(row['GERP++_RS'])
            scores.append(min(1.0, max(0, gerp / 6)))
        
        return np.mean(scores) if scores else 0.5

    def score_splicing(self, row: pd.Series) -> float:
        """Score variant based on splicing predictions"""
        scores = []
        
        if pd.notna(row['dbscSNV_ADA_SCORE']) and row['dbscSNV_ADA_SCORE'] != '.':
            ada = float(row['dbscSNV_ADA_SCORE'])
            scores.append(1.0 if ada >= 0.6 else 0.0)
        
        if pd.notna(row['dbscSNV_RF_SCORE']) and row['dbscSNV_RF_SCORE'] != '.':
            rf = float(row['dbscSNV_RF_SCORE'])
            scores.append(1.0 if rf >= 0.6 else 0.0)
        
        return np.mean(scores) if scores else 0.5
        
    def get_dynamic_weights(self, row: pd.Series) -> Dict[str, float]:
        """Dynamically adjust weights based on variant characteristics"""
        # Start with default weights
        weights = self.weights.copy()
        
        # Check for important variant characteristics
        has_pvs1 = self._has_pvs1_evidence(row)
        is_high_impact = self._is_high_impact_variant(row)
        in_high_impact_gene = self._is_in_high_impact_gene(row)
        
        # Dynamic adjustment for PVS1=1 variants
        if has_pvs1 and is_high_impact:
            if in_high_impact_gene:
                # For critical variants in high-impact genes, increase clinical and functional importance
                weights['clinical_evidence'] = 0.45   # Increase from 0.40
                weights['functional_impact'] = 0.20   # Increase from 0.15
                weights['population_freq'] = 0.10     # Decrease from 0.15 (less penalty for frequency)
            else:
                # For other PVS1 variants, still adjust but less dramatically
                weights['clinical_evidence'] = 0.42   # Slight increase
                weights['functional_impact'] = 0.18   # Slight increase
                weights['population_freq'] = 0.12     # Slight decrease
        
        # Normalize weights to ensure they sum to 1
        total = sum(weights.values())
        weights = {k: v/total for k, v in weights.items()}
        
        return weights

    def calculate_final_score(self, row: pd.Series) -> Dict[str, float]:
        """Calculate final weighted score with component scores and minimum thresholds"""
        # Calculate component scores - pass row to frequency scoring for context awareness
        component_scores = {
            'functional_impact': self.score_functional_impact(row),
            'population_freq': self.score_population_frequency(row[self.column_names['freq_gnomad']], row),
            'clinical_evidence': self.score_clinical_evidence(row),
            'computational': self.score_computational_predictions(row),
            'conservation': self.score_conservation(row),
            'splicing': self.score_splicing(row)
        }
        
        # Get dynamically adjusted weights based on variant characteristics
        weights = self.get_dynamic_weights(row)
        
        # Calculate weighted sum
        final_score = sum(score * weights[component] 
                         for component, score in component_scores.items())
        
        # Apply minimum score thresholds for critical variants
        has_pvs1 = self._has_pvs1_evidence(row)
        is_high_impact = self._is_high_impact_variant(row)
        in_high_impact_gene = self._is_in_high_impact_gene(row)
        
        # Set minimum score floor based on variant criticality
        if has_pvs1 and is_high_impact:
            if in_high_impact_gene:
                # Critical gene PVS1 variants get at least Tier 2 ranking
                min_score = 0.65  # High end of Tier 2
                final_score = max(final_score, min_score)
            else:
                # Other gene PVS1 variants get at least low Tier 2 ranking
                min_score = 0.60  # Low end of Tier 2
                final_score = max(final_score, min_score)
        
        # Add final score to component scores
        component_scores['final_score'] = final_score
        
        return component_scores

    def classify_variant(self, score: float) -> str:
        """Refined classification thresholds"""
        if score >= 0.75:            # Lowered threshold for Tier 1
            return "Tier 1 - High Priority"
        elif score >= 0.60:          # Medium priority threshold
            return "Tier 2 - Medium Priority"
        elif score >= 0.45:          # Lower threshold for Tier 3
            return "Tier 3 - Low Priority"
        else:
            return "Tier 4 - Benign/Likely Benign"

    def process_variants(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process all variants in the dataframe"""
        # Create new columns for scores and classification
        df['component_scores'] = df.apply(self.calculate_final_score, axis=1)
        df['final_score'] = df['component_scores'].apply(lambda x: x['final_score'])
        df['classification'] = df['final_score'].apply(self.classify_variant)
        
        # Sort by final score
        df_sorted = df.sort_values('final_score', ascending=False)
        
        return df_sorted
