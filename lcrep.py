import pandas as pd
import jinja2
import os
import sys
from datetime import datetime
import re

class ClinicalReportGenerator:
    def __init__(self):
        self.template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Variant Analysis Report</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/3.7.0/chart.min.js"></script>
    <style>
        :root {
            --primary-color: #2c3e50;
            --secondary-color: #34495e;
            --accent-color: #3498db;
            --success-color: #27ae60;
            --warning-color: #f1c40f;
            --danger-color: #e74c3c;
            --light-gray: #f5f6fa;
            --alternate-row: #f8fafc;
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        }
        
        body {
            background-color: var(--light-gray);
            color: var(--primary-color);
            line-height: 1.6;
        }
        
        .container {
            max-width: 1800px;
            margin: 0 auto;
            padding: 20px;
        }
        
        .header {
            background-color: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        
        .header-content {
            flex: 1;
        }
        
        .header h1 {
            color: var(--primary-color);
            margin-bottom: 10px;
        }

        .logo {
            text-align: right;
        }

        .logo img {
            height: 64px;
            width: auto;
        }

        .database-dropdown {
            position: relative;
            display: inline-block;
            padding: 4px 8px;
        }

        .database-dropdown::before {
            content: '';
            position: absolute;
            top: -30px;
            right: -30px;
            bottom: -30px;
            left: -200px;
            z-index: -1;
        }

        .database-btn {
            color: var(--accent-color);
            text-decoration: none;
            font-weight: bold;
            cursor: pointer;
            font-size: 15px;
            padding: 4px 8px;
            position: relative;
        }

        .database-dropdown::after {
            content: '';
            position: absolute;
            top: 0;
            right: 100%;
            height: 100%;
            width: 20px;
            z-index: 999;
        }

        .dropdown-content {
            visibility: hidden;
            opacity: 0;
            position: absolute;
            background-color: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.2);
            border-radius: 4px;
            right: 100%;
            top: 50%;
            transform: translateY(-50%);
            margin-right: 8px;
            white-space: nowrap;
            padding: 8px 16px;
            z-index: 1000;
            transition: visibility 0s, opacity 0.2s ease-in-out;
        }

        .database-dropdown:hover .dropdown-content,
        .dropdown-content:hover {
            visibility: visible;
            opacity: 1;
        }

        .dropdown-content::after {
            content: '';
            position: absolute;
            top: -20px;
            right: -20px;
            bottom: -20px;
            left: -20px;
            z-index: -1;
        }

        .dropdown-content {
            display: flex;
            gap: 16px;
        }

        .dropdown-content a {
            color: var(--accent-color);
            text-decoration: none;
            font-size: 15px;
            font-weight: 500;
            padding: 4px 0;
            white-space: nowrap;
        }

        .dropdown-content a:hover {
            text-decoration: underline;
        }

        .dropdown-content a:not(:last-child):after {
            content: '|';
            color: #ccc;
            margin-left: 16px;
            font-size: 15px;
        }

        .controls {
            background-color: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }

        .search-box {
            width: 100%;
            padding: 8px;
            margin-bottom: 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }

        .column-toggle {
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            margin-bottom: 10px;
        }

        .column-toggle label {
            display: flex;
            align-items: center;
            gap: 5px;
            font-size: 14px;
        }

        .stats-chart-container {
            background-color: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            display: grid;
            grid-template-columns: 3fr 2fr;
            gap: 20px;
        }

        .metrics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
        }

        .metric-card {
            background-color: var(--light-gray);
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }

        .metric-value {
            font-size: 1.5em;
            font-weight: bold;
            color: var(--primary-color);
        }

        .metric-label {
            color: var(--secondary-color);
            font-size: 0.9em;
            margin-top: 5px;
        }

        .chart-container {
            height: 100%;
            min-height: 300px;
            display: flex;
            justify-content: center;
            align-items: center;
        }

        .variants-table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9em;
            table-layout: fixed;
        }

        .col-variant { width: 140px; }
        .col-gene { width: 70px; }
        .col-dbsnp { width: 90px; }
        .col-cosmic { width: 90px; }
        .col-location { width: 100px; }
        .col-effect { width: 120px; }
        .col-clnhgvs { width: 110px; }
        .col-quality { width: 70px; }
        .col-significance { width: 150px; }
        .col-acmg-class { width: 150px; }
        .col-acmg-criteria { width: 150px; }
        .col-database { width: 90px; }
        .col-gnomad-af { width: 80px; text-align: right; }
        .col-depth { width: 80px; text-align: right; }
        .col-disease { width: 120px; }

        .not-available {
            color: #666;
            font-style: italic;
        }

        .variants-table th,
        .variants-table td {
            padding: 12px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            text-align: left;
        }

        .variants-table td:hover {
            overflow: visible;
            white-space: normal;
            word-break: break-word;
            position: relative;
            background-color: white;
            z-index: 1;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .variants-table th {
            background-color: var(--primary-color);
            color: white;
            font-weight: 600;
            cursor: pointer;
            position: sticky;
            top: 0;
            text-align: left;
        }

        .variants-table th:hover {
            background-color: var(--secondary-color);
        }

        .variants-table td {
            border-bottom: 1px solid #eee;
        }

        .variants-table tr:nth-child(even) {
            background-color: var(--alternate-row);
        }

        .variants-table tr:hover {
            background-color: #e8f4fc;
        }

        .significance-badge {
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: 600;
            display: inline-block;
        }

        .pathogenic {
            background-color: #FF0000;
            color: white;
        }

        .likely-pathogenic {
            background-color: #FF6666;
            color: white;
        }

        .uncertain {
            background-color: #ADD8E6;
            color: #2c3e50;
        }

        .likely-benign {
            background-color: #4CAF50;
            color: white;
        }

        .benign {
            background-color: #4CAF50;
            color: white;
        }

        .conflicting {
            background-color: #FFA500;
            color: white;
        }

        .acmg-badge {
            padding: 4px 8px;
            border-radius: 4px;
            font-weight: 600;
            display: inline-block;
        }

        .acmg-pathogenic {
            background-color: #FF0000;
            color: white;
        }

        .acmg-likely-pathogenic {
            background-color: #FF6666;
            color: white;
        }

        .acmg-uncertain {
            background-color: #BFDFED;
            color: #2c3e50;
        }

        .acmg-likely-benign {
            background-color: #66BB6A;
            color: white;
        }

        .acmg-benign {
            background-color: #81C784;
            color: white;
        }

        .acmg-criteria {
            font-family: monospace;
            font-size: 1.2em;
            color: #666;
            white-space: pre-wrap;
            font-weight: bold;
        }

        .filter-pass {
            color: var(--success-color);
            font-weight: bold;
        }

        .filter-fail {
            color: var(--danger-color);
            font-weight: bold;
        }

        .gene-link {
            color: var(--accent-color);
            text-decoration: none;
            font-weight: bold;
        }

        .gene-link:hover {
            text-decoration: underline;
        }

        .footer {
            background-color: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-top: 20px;
            text-align: center;
            color: var(--secondary-color);
            font-size: 0.9em;
        }

        @media (max-width: 768px) {
            .container {
                padding: 10px;
            }

            .variants-table {
                display: block;
                overflow-x: auto;
            }

            .header {
                flex-direction: column;
                text-align: center;
            }

            .logo {
                margin-top: 10px;
                text-align: center;
            }
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="header-content">
                <h1>Whole Exome Sequencing Analysis Report - Sample ID: {{ sample_id }}</h1>
                <p>Generated on {{ generation_date }}</p>
            </div>
            <div class="logo">
                <img src="lifecode.png" alt="LifeCode Logo" style="height: 64px; width: auto;">
            </div>
        </div>

        <div class="controls">
            <input type="text" class="search-box" id="variantSearch" placeholder="Search for genes, variants, or diseases...">
            <div class="column-toggle" id="columnToggle">
            </div>
        </div>

        <div class="stats-chart-container">
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-value">{{ total_variants }}</div>
                    <div class="metric-label">Total Variants</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ pathogenic_count }}</div>
                    <div class="metric-label">Pathogenic/Likely Pathogenic</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ vus_count }}</div>
                    <div class="metric-label">Variants of Uncertain Significance</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ benign_count }}</div>
                    <div class="metric-label">Benign/Likely Benign</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{{ pass_count }}</div>
                    <div class="metric-label">PASS Filter</div>
                </div>
            </div>
            <div class="chart-container">
                <canvas id="variantDistribution"></canvas>
            </div>
        </div>

        <div class="variants-container">
                <table class="variants-table" id="variantsTable">
                    <thead>
                        <tr>
                            <th class="col-variant" data-sort="variant">Variant</th>
                            <th class="col-gene" data-sort="gene">Gene</th>
                            <th class="col-dbsnp" data-sort="dbsnp">dbSNP ID</th>
                            <th class="col-cosmic" data-sort="cosmic">COSMIC</th>
                            <th class="col-location" data-sort="location">Location</th>
                            <th class="col-effect" data-sort="effect">Effect</th>
                            <th class="col-clnhgvs" data-sort="clnhgvs">CLNHGVS</th>
                            <th class="col-quality" data-sort="quality">Quality</th>
                            <th class="col-significance" data-sort="significance">Clinical Significance</th>
                            <th class="col-acmg-class" data-sort="acmg">ACMG Classification</th>
                            <th class="col-acmg-criteria" data-sort="criteria">ACMG Met Criteria</th>
                            <th class="col-database" data-sort="database">Database</th>
                            <th class="col-gnomad-af" data-sort="gnomad">gnomAD AF</th>
                            <th class="col-disease" data-sort="disease">Disease</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for variant in variants %}
                        <tr>
                            <td class="col-variant">{{ variant.Variant }}</td>
                            <td class="col-gene">
                                <a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene={{ variant.Gene }}" 
                                   class="gene-link" target="_blank" 
                                   title="View {{ variant.Gene }} in GeneCards">{{ variant.Gene }}</a>
                            </td>
                            <td class="col-dbsnp">{{ variant.dbSNP_ID }}</td>
                            <td class="col-cosmic">{{ variant.COSMIC_ID }}</td>
                            <td class="col-location">{{ variant.Location }}</td>
                            <td class="col-effect">{{ variant.Effect }}</td>
                            <td class="col-clnhgvs">{{ variant.CLNHGVS }}</td>
                            <td class="col-quality" title="{{ variant.Quality }}">
                                <span class="{{ variant.quality_class }}">
                                    {{ variant.Quality }}
                                </span>
                            </td>
                            <td class="col-significance">
                                <span class="significance-badge {{ variant.significance_class }}">
                                    {{ variant.Clinical_Significance }}
                                </span>
                            </td>
                            <td class="col-acmg-class">
                                {% if variant.ACMG_Classification %}
                                <span class="acmg-badge {{ variant.acmg_class }}">
                                    {{ variant.ACMG_Classification }}
                                </span>
                                {% else %}
                                <span class="not-available">Not Available</span>
                                {% endif %}
                            </td>
                            <td class="col-acmg-criteria">
                                <span class="acmg-criteria">{{ variant.ACMG_Met_Criteria }}</span>
                            </td>
                            <td class="col-database">
                                <div class="database-dropdown">
                                    <span class="database-btn">View</span>
                                    <div class="dropdown-content">
                                        {% for db_name, db_url in variant.database_links.items() %}
                                        <a href="{{ db_url }}" target="_blank">{{ db_name }}</a>
                                        {% endfor %}
                                    </div>
                                </div>
                            </td>
                            <td class="col-gnomad-af">
                                {% if variant.gnomAD_AF_formatted != "N/A" %}
                                    {{ variant.gnomAD_AF_formatted }}
                                {% else %}
                                    <span class="not-available">Not Available</span>
                                {% endif %}
                            </td>
                            <td class="col-disease">{{ variant.Disease }}</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
            
            <div class="footer">
                <p>This report contains the top 300 prioritized variants. For full results, please refer to the complete data file.</p>
                <p><small>Tip: Sort columns by clicking headers. Priority sorting implemented for pathogenicity.</small></p>
            </div>
        </div>

        <script>
            // Initialize variant distribution chart
            const ctx = document.getElementById('variantDistribution').getContext('2d');
            new Chart(ctx, {
                type: 'doughnut',
                data: {
                    labels: ['Pathogenic/Likely Pathogenic', 'VUS', 'Benign/Likely Benign', 'Conflicting'],
                    datasets: [{
                        data: [{{ pathogenic_count }}, {{ vus_count }}, {{ benign_count }}, {{ conflicting_count }}],
                        backgroundColor: ['#FF0000', '#ADD8E6', '#4CAF50', '#FFA500']
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {
                        legend: {
                            position: 'right',
                            labels: {
                                boxWidth: 12,
                                padding: 15
                            }
                        },
                        title: {
                            display: true,
                            text: 'Variant Classification Distribution',
                            padding: {
                                top: 10,
                                bottom: 20
                            }
                        }
                    }
                }
            });

            // Column visibility toggle
            const columnToggle = document.getElementById('columnToggle');
            const table = document.getElementById('variantsTable');
            const columns = [
                {name: 'Variant', class: 'col-variant'},
                {name: 'Gene', class: 'col-gene'},
                {name: 'dbSNP ID', class: 'col-dbsnp'},
                {name: 'COSMIC', class: 'col-cosmic'},
                {name: 'Location', class: 'col-location'},
                {name: 'Effect', class: 'col-effect'},
                {name: 'CLNHGVS', class: 'col-clnhgvs'},
                {name: 'Quality', class: 'col-quality'},
                {name: 'Clinical Significance', class: 'col-significance'},
                {name: 'ACMG Classification', class: 'col-acmg-class'},
                {name: 'ACMG Met Criteria', class: 'col-acmg-criteria'},
                {name: 'Database', class: 'col-database'},
                {name: 'gnomAD AF', class: 'col-gnomad-af'},
                {name: 'Disease', class: 'col-disease'}
            ];

            columns.forEach(col => {
                const label = document.createElement('label');
                const checkbox = document.createElement('input');
                checkbox.type = 'checkbox';
                checkbox.checked = true;
                checkbox.addEventListener('change', () => {
                    const cells = table.getElementsByClassName(col.class);
                    for (let cell of cells) {
                        cell.style.display = checkbox.checked ? '' : 'none';
                    }
                });
                label.appendChild(checkbox);
                label.appendChild(document.createTextNode(col.name));
                columnToggle.appendChild(label);
            });

            // Search functionality
            const searchBox = document.getElementById('variantSearch');
            searchBox.addEventListener('input', () => {
                const searchTerm = searchBox.value.toLowerCase();
                const rows = table.getElementsByTagName('tr');

                for (let i = 1; i < rows.length; i++) {
                    const row = rows[i];
                    const text = row.textContent.toLowerCase();
                    row.style.display = text.includes(searchTerm) ? '' : 'none';
                }
            });

            // Add significance and ACMG priority mapping
            const significancePriority = {
                'pathogenic': 1,
                'likely-pathogenic': 2,
                'conflicting': 3,
                'uncertain': 4,
                'likely-benign': 5,
                'benign': 6,
                'unknown': 7
            };

            const acmgPriority = {
                'acmg-pathogenic': 1,
                'acmg-likely-pathogenic': 2,
                'acmg-uncertain': 3,
                'acmg-likely-benign': 4,
                'acmg-benign': 5,
                'unknown': 6
            };

            // Function to get significance priority
            function getSignificancePriority(element) {
                const badge = element.querySelector('.significance-badge');
                if (!badge) return 999;
                
                const classes = Array.from(badge.classList);
                const significanceClass = classes.find(cls => significancePriority.hasOwnProperty(cls));
                
                return significanceClass ? significancePriority[significanceClass] : 999;
            }

            // Function to get ACMG priority
            function getACMGPriority(element) {
                const badge = element.querySelector('.acmg-badge');
                if (!badge) return 999;
                
                const classes = Array.from(badge.classList);
                const acmgClass = classes.find(cls => acmgPriority.hasOwnProperty(cls));
                
                return acmgClass ? acmgPriority[acmgClass] : 999;
            }

            // Function to parse numeric values including scientific notation
            function parseNumericValue(value) {
                if (value === null || value === undefined || value === '') return NaN;
                if (typeof value === 'number') return value;
                if (value.toLowerCase() === 'not available') return NaN;
                return parseFloat(value);
            }

            // Sorting functionality
            let currentSort = { column: null, ascending: true };

            function sortTable(column) {
                const tbody = table.getElementsByTagName('tbody')[0];
                const rows = Array.from(tbody.getElementsByTagName('tr'));

                if (currentSort.column === column) {
                    currentSort.ascending = !currentSort.ascending;
                } else {
                    currentSort = { column: column, ascending: true };
                }

                rows.sort((a, b) => {
                    if (column === 'significance') {
                        const aVal = getSignificancePriority(a.querySelector('.col-significance'));
                        const bVal = getSignificancePriority(b.querySelector('.col-significance'));
                        return currentSort.ascending ? aVal - bVal : bVal - aVal;
                    }

                    if (column === 'acmg') {
                        const aVal = getACMGPriority(a.querySelector('.col-acmg-class'));
                        const bVal = getACMGPriority(b.querySelector('.col-acmg-class'));
                        return currentSort.ascending ? aVal - bVal : bVal - aVal;
                    }
                    
                    if (column === 'gnomad') {
                        const aText = a.querySelector(`.col-${column}-af`).textContent;
                        const bText = b.querySelector(`.col-${column}-af`).textContent;
                        const aVal = parseNumericValue(aText);
                        const bVal = parseNumericValue(bText);
                        
                        // Handle cases where one or both values are NaN
                        if (isNaN(aVal) && isNaN(bVal)) return 0;
                        if (isNaN(aVal)) return currentSort.ascending ? 1 : -1;
                        if (isNaN(bVal)) return currentSort.ascending ? -1 : 1;
                        
                        return currentSort.ascending ? aVal - bVal : bVal - aVal;
                    }
                    
                    const aVal = a.querySelector(`.col-${column}`).textContent;
                    const bVal = b.querySelector(`.col-${column}`).textContent;
                    return currentSort.ascending ? 
                        aVal.localeCompare(bVal, undefined, {numeric: true}) : 
                        bVal.localeCompare(aVal, undefined, {numeric: true});
                });

                rows.forEach(row => tbody.removeChild(row));
                rows.forEach(row => tbody.appendChild(row));
            }

            // Add click handlers for sorting
            const headers = table.getElementsByTagName('th');
            for (let header of headers) {
                const sortColumn = header.getAttribute('data-sort');
                if (sortColumn) {
                    header.addEventListener('click', () => sortTable(sortColumn));
                }
            }
        </script>
    </body>
</html>"""

    def get_significance_class(self, significance):
        """Determine CSS class based on clinical significance."""
        if pd.isna(significance):
            return 'unknown'
        
        significance = str(significance).lower()
        if 'conflicting' in significance:
            return 'conflicting'
        elif 'pathogenic' in significance and 'likely' not in significance:
            return 'pathogenic'
        elif 'likely_pathogenic' in significance or 'likely pathogenic' in significance:
            return 'likely-pathogenic'
        elif 'uncertain_significance' in significance or 'uncertain significance' in significance:
            return 'uncertain'
        elif 'likely_benign' in significance or 'likely benign' in significance:
            return 'likely-benign'
        elif 'benign' in significance:
            return 'benign'
        return 'unknown'

    def get_acmg_class(self, classification):
        """Determine CSS class based on ACMG classification."""
        if pd.isna(classification):
            return 'unknown'
        
        classification = str(classification).lower()
        if 'pathogenic' in classification and 'likely' not in classification:
            return 'acmg-pathogenic'
        elif 'likely pathogenic' in classification:
            return 'acmg-likely-pathogenic'
        elif 'uncertain' in classification:
            return 'acmg-uncertain'
        elif 'likely benign' in classification:
            return 'acmg-likely-benign'
        elif 'benign' in classification:
            return 'acmg-benign'
        return 'unknown'

    def get_quality_class(self, quality):
        """Determine CSS class based on quality value."""
        if pd.isna(quality):
            return 'filter-fail'
        return 'filter-pass' if quality == 'PASS' else 'filter-fail'

    def format_af(self, af_value):
        """Format allele frequency value."""
        try:
            if pd.isna(af_value) or af_value == '.':
                return "N/A"
            
            # Handle scientific notation
            if isinstance(af_value, str) and 'e' in af_value.lower():
                af = float(af_value)
            else:
                af = float(af_value)
                
            if af == 0:
                return "0.0000"
            return f"{af:.4f}"
        except (ValueError, TypeError):
            return "N/A"

    def get_database_links(self, gene, variant_dict):
        """Generate database links for a given gene and variant."""
        # Extract necessary components from variant information
        chrom = variant_dict.get('Chr', '').replace('chr', '')
        pos = str(variant_dict.get('Start', ''))
        ref = variant_dict.get('Ref', '')
        alt = variant_dict.get('Alt', '')
        variant_str = f"chr{chrom}-{pos}-{ref}-{alt}-hg19"
        hgvs_g = f"{variant_dict.get('Chr')}:g.{pos}{ref}>{alt}"
        
        # Determine ClinVar term priority: CLNHGVS > dbSNP ID
        clinvar_term = None
        if variant_dict.get('CLNHGVS') and not pd.isna(variant_dict.get('CLNHGVS')):
            clinvar_term = variant_dict['CLNHGVS']
        elif variant_dict.get('rs') and not pd.isna(variant_dict.get('rs')) and variant_dict['rs'] != '.':
            clinvar_term = variant_dict['rs']
        
        # Create database links
        links = {
            'ClinVar': f"https://www.ncbi.nlm.nih.gov/clinvar/?term={clinvar_term}" if clinvar_term else "#",
            'Franklin': f"https://franklin.genoox.com/clinical-db/variant/snp/{variant_str}",
            'Varsome': f"https://varsome.com/variant/hg19/{variant_dict.get('rs') if variant_dict.get('rs') and variant_dict.get('rs') != '.' else variant_str}?annotation-mode=germline",
            'OncoKB': f"https://www.oncokb.org/hgvsg/{hgvs_g}?refGenome=GRCh37"
        }
        
        # Add COSMIC link if LEGACY_ID exists and contains COSM
        if variant_dict.get('LEGACY_ID') and not pd.isna(variant_dict.get('LEGACY_ID')):
            legacy_id = str(variant_dict['LEGACY_ID'])
            if 'COSM' in legacy_id:
                cosmic_id = legacy_id.replace('COSM', '')  # Remove 'COSM' prefix if present
                links['COSMIC'] = f"https://cancer.sanger.ac.uk/cosmic/mutation/overview?id={cosmic_id}"
        
        # Add OMIM link if disease info contains an OMIM number
        if variant_dict.get('CLNDN') and not pd.isna(variant_dict.get('CLNDN')):
            links['OMIM'] = f"https://omim.org/search/?search={variant_dict['CLNDN']}"
            
        # Add Orphanet link if we have Orphanet data
        if variant_dict.get('Orpha') and not pd.isna(variant_dict.get('Orpha')):
            # Extract orphanet ID from Orpha field if present
            orpha_match = re.search(r'(\d+)\|', str(variant_dict['Orpha']))
            if orpha_match:
                orpha_id = orpha_match.group(1)
                links['Orphanet'] = f"https://www.orpha.net/consor/cgi-bin/OC_Exp.php?Lng=EN&Expert={orpha_id}"
        
        return links

    def count_variants_by_significance(self, df):
        """Count variants by their clinical significance."""
        pathogenic_count = 0
        benign_count = 0
        vus_count = 0
        conflicting_count = 0

        for idx, row in df.iterrows():
            # Check CLNSIG first, then fall back to ACMG_Classification if needed
            significance = row['CLNSIG']
            if pd.isna(significance) or significance == '':
                significance = row['ACMG_Classification']
                if pd.isna(significance) or significance == '':
                    continue
            
            significance = str(significance).lower()
            if 'conflicting' in significance:
                conflicting_count += 1
            elif ('pathogenic' in significance and 'likely' not in significance) or \
                 'likely_pathogenic' in significance or 'likely pathogenic' in significance:
                pathogenic_count += 1
            elif 'uncertain_significance' in significance or 'uncertain significance' in significance:
                vus_count += 1
            elif 'benign' in significance:
                benign_count += 1

        return pathogenic_count, benign_count, vus_count, conflicting_count

    def generate_vcf(self, df, output_vcf):
        """Generate VCF file from the dataframe."""
        # Write VCF header
        with open(output_vcf, 'w') as f:
            # Write VCF metadata
            f.write('##fileformat=VCFv4.2\n')
            f.write(f'##fileDate={datetime.now().strftime("%Y%m%d")}\n')
            f.write('##source=BioWESReport\n')
            f.write('##reference=hg19\n')
            
            # Write format headers for additional fields
            f.write('##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="Population Allele Frequency (gnomAD)">\n')
            f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
            f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Symbol">\n')
            f.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical Significance">\n')
            f.write('##INFO=<ID=ACMG,Number=.,Type=String,Description="ACMG Classification">\n')
            f.write('##INFO=<ID=EFFECT,Number=1,Type=String,Description="Variant Effect">\n')
            f.write('##INFO=<ID=LOC,Number=1,Type=String,Description="Variant Location">\n')
            
            # Write VCF column headers
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            
            # Write variant data
            for _, row in df.iterrows():
                # Prepare INFO field
                info_fields = []
                
                if not pd.isna(row['Freq_gnomAD_genome_ALL']) and row['Freq_gnomAD_genome_ALL'] != '.':
                    info_fields.append(f"gnomAD_AF={row['Freq_gnomAD_genome_ALL']}")
                
                if not pd.isna(row['DP']):
                    try:
                        dp = int(row['DP'])
                        info_fields.append(f"DP={dp}")
                    except (ValueError, TypeError):
                        pass
                
                if not pd.isna(row['Ref.Gene']):
                    info_fields.append(f"GENE={row['Ref.Gene']}")
                
                if not pd.isna(row['CLNSIG']):
                    info_fields.append(f"CLNSIG={row['CLNSIG']}")
                
                if not pd.isna(row['ACMG_Classification']):
                    info_fields.append(f"ACMG={row['ACMG_Classification']}")
                
                if not pd.isna(row['effect']):
                    info_fields.append(f"EFFECT={row['effect']}")
                
                if not pd.isna(row['location']):
                    info_fields.append(f"LOC={row['location']}")
                
                info = ';'.join(info_fields)
                
                # Write variant line
                chrom = str(row['Chr']).replace('chr', '')  # Remove 'chr' prefix if present
                
                variant_line = [
                    chrom,
                    str(row['Start']),
                    str(row['rs']) if not pd.isna(row['rs']) and row['rs'] != '.' else '.',
                    str(row['Ref']),
                    str(row['Alt']),
                    '.',  # QUAL
                    str(row['quality']) if not pd.isna(row['quality']) else '.',
                    info
                ]
                f.write('\t'.join(variant_line) + '\n')

    def generate_report(self, input_file, output_file, sample_id=None):
        """Generate HTML report and VCF from TSV/Excel file."""
        # Determine file type and read accordingly
        file_ext = os.path.splitext(input_file)[1].lower()
        
        if file_ext == '.xlsx' or file_ext == '.xls':
            df = pd.read_excel(input_file)
        elif file_ext == '.tsv' or file_ext == '.txt':
            df = pd.read_csv(input_file, sep='\t')
        elif file_ext == '.csv':
            df = pd.read_csv(input_file)
        else:
            raise ValueError(f"Unsupported file format: {file_ext}")
        
        # Extract sample ID from filename if not provided
        if not sample_id:
            sample_id = os.path.splitext(os.path.basename(input_file))[0]
        
        # Ensure all required columns exist
        required_columns = ['Chr', 'Start', 'Ref', 'Alt', 'Ref.Gene', 'location', 
                          'rs', 'LEGACY_ID', 'quality', 'CLNSIG', 'ACMG_Classification', 
                          'ACMG_Met_Criteria', 'CLNHGVS', 'CLNDN', 'Freq_gnomAD_genome_ALL', 
                          'effect', 'DP']
        
        for col in required_columns:
            if col not in df.columns:
                df[col] = None  # Add empty column if missing
        
        # Limit to top 300 variants if there are more
        df = df.head(300)
        
        # Generate VCF file
        output_vcf = output_file.rsplit('.', 1)[0] + '.vcf'
        self.generate_vcf(df, output_vcf)
        
        # Count variants by significance
        pathogenic_count, benign_count, vus_count, conflicting_count = self.count_variants_by_significance(df)
        
        # Add CSS classes and database links for styling
        variants = []
        for _, row in df.iterrows():
            variant_dict = row.to_dict()
            
            # Create formatted variant string
            variant_dict['Variant'] = f"{row['Chr']}:{row['Start']} {row['Ref']}->{row['Alt']}"
            
            # Map columns to their report equivalents
            variant_dict['Gene'] = variant_dict['Ref.Gene']
            variant_dict['dbSNP_ID'] = variant_dict['rs'] if not pd.isna(variant_dict['rs']) and variant_dict['rs'] != '.' else 'Novel'
            
            # Extract COSMIC ID from LEGACY_ID if it contains "COSM"
            cosmic_id = "Not Available"
            if not pd.isna(variant_dict['LEGACY_ID']) and 'COSM' in str(variant_dict['LEGACY_ID']):
                cosmic_id = variant_dict['LEGACY_ID']
            variant_dict['COSMIC_ID'] = cosmic_id
            
            variant_dict['Location'] = variant_dict['location']
            variant_dict['Effect'] = variant_dict['effect']
            variant_dict['Quality'] = variant_dict['quality']
            variant_dict['Disease'] = variant_dict['CLNDN'] if not pd.isna(variant_dict['CLNDN']) else ''
            
            # Clinical significance - prefer CLNSIG but fallback to ACMG_Classification if needed
            clinical_sig = variant_dict['CLNSIG']
            if pd.isna(clinical_sig) or clinical_sig == '':
                clinical_sig = variant_dict['ACMG_Classification']
            variant_dict['Clinical_Significance'] = clinical_sig if not pd.isna(clinical_sig) else 'Not Available'
            
            # Add classes and formatting
            variant_dict['significance_class'] = self.get_significance_class(clinical_sig)
            variant_dict['acmg_class'] = self.get_acmg_class(variant_dict['ACMG_Classification'])
            variant_dict['quality_class'] = self.get_quality_class(variant_dict['quality'])
            variant_dict['database_links'] = self.get_database_links(variant_dict['Gene'], variant_dict)
            
            # Format gnomAD AF value
            variant_dict['gnomAD_AF_formatted'] = self.format_af(variant_dict['Freq_gnomAD_genome_ALL'])
            
            # Convert depth to int if possible
            try:
                variant_dict['Depth'] = int(variant_dict['DP']) if not pd.isna(variant_dict['DP']) else 0
            except (ValueError, TypeError):
                variant_dict['Depth'] = 0
            
            variants.append(variant_dict)
        
        # Prepare template data
        template_data = {
            'sample_id': sample_id,
            'generation_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'total_variants': len(df),
            'pathogenic_count': pathogenic_count,
            'vus_count': vus_count,
            'benign_count': benign_count,
            'conflicting_count': conflicting_count,
            'pass_count': df['quality'].str.contains('PASS', case=True, na=False).sum(),
            'variants': variants
        }
        
        # Render template
        template = jinja2.Template(self.template)
        html_content = template.render(**template_data)
        
        # Write HTML file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"HTML report generated successfully: {output_file}")
        print(f"VCF file generated successfully: {output_vcf}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python bio_WES_report_v0.2.py input.tsv output.html [sample_id]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sample_id = sys.argv[3] if len(sys.argv) > 3 else None
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist")
        sys.exit(1)
    
    if not output_file.lower().endswith('.html'):
        output_file += '.html'
    
    try:
        generator = ClinicalReportGenerator()
        generator.generate_report(input_file, output_file, sample_id)
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
