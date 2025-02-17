import pandas as pd
import jinja2
import os
import sys
from datetime import datetime
import json

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

        .col-chr { width: 55px; }
        .col-pos { width: 80px; }
        .col-ref { width: 45px; }
        .col-alt { width: 45px; }
        .col-gene { width: 70px; }
        .col-location { width: 100px; }
        .col-dbsnp { width: 90px; }
        .col-cosmic { width: 90px; }
        .col-transcript { width: 100px; }
        .col-clnhgvs { width: 100px; }
        .col-filter { width: 70px; }
        .col-significance { width: 150px; }
        .col-acmg-class { width: 150px; }
        .col-acmg-criteria { width: 150px; }
        .col-database { width: 90px; }
        .col-af { width: 60px; text-align: right; }
        .col-gnomad-af { width: 80px; text-align: right; }
        .col-depth { width: 80px; text-align: right; }
        .col-disease { width: 100px; }

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
                            <th class="col-chr" data-sort="chr">Chr</th>
                            <th class="col-pos" data-sort="pos">Position</th>
                            <th class="col-ref" data-sort="ref">Ref</th>
                            <th class="col-alt" data-sort="alt">Alt</th>
                            <th class="col-gene" data-sort="gene">Gene</th>
                            <th class="col-location" data-sort="location">Location</th>
                            <th class="col-dbsnp" data-sort="dbsnp">dbSNP ID</th>
                            <th class="col-cosmic" data-sort="cosmic">COSMIC</th>
                            <th class="col-transcript" data-sort="transcript">Transcript</th>
                            <th class="col-clnhgvs" data-sort="clnhgvs">CLNHGVS</th>
                            <th class="col-filter" data-sort="filter">Filter</th>
                            <th class="col-significance" data-sort="significance">Clinical Significance</th>
                            <th class="col-acmg-class" data-sort="acmg">ACMG Classification</th>
                            <th class="col-acmg-criteria" data-sort="criteria">ACMG Met Criteria</th>
                            <th class="col-database" data-sort="database">Database</th>
                            <th class="col-af" data-sort="af">AF</th>
                            <th class="col-gnomad-af" data-sort="gnomad">P-AF</th>
                            <th class="col-depth" data-sort="depth">Depth</th>
                            <th class="col-disease" data-sort="disease">Disease</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for variant in variants %}
                        <tr>
                            <td class="col-chr">{{ variant.CHROM }}</td>
                            <td class="col-pos">{{ variant.POS }}</td>
                            <td class="col-ref">{{ variant.REF }}</td>
                            <td class="col-alt">{{ variant.ALT }}</td>
                            <td class="col-gene">
                                <a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene={{ variant.Gene }}" 
                                   class="gene-link" target="_blank" 
                                   title="View {{ variant.Gene }} in GeneCards">{{ variant.Gene }}</a>
                            </td>
                            <td class="col-location">{{ variant.Location }}</td>
                            <td class="col-dbsnp">{{ variant.dbSNP_ID }}</td>
                            <td class="col-cosmic">{{ variant.COSMIC_ID }}</td>
                            <td class="col-transcript">{{ variant.Transcript }}</td>
                            <td class="col-clnhgvs">{{ variant.CLNHGVS }}</td>
                            <td class="col-filter" title="{{ variant.FILTER }}">
                                <span class="{{ variant.filter_class }}">
                                    {{ variant.FILTER }}
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
                            <td class="col-af">
                                {% if variant.AF_formatted != "N/A" %}
                                    {{ variant.AF_formatted }}
                                {% else %}
                                    <span class="not-available">Not Available</span>
                                {% endif %}
                            </td>
                            <td class="col-gnomad-af">
                                {% if variant.gnomAD_AF_formatted != "N/A" %}
                                    {{ variant.gnomAD_AF_formatted }}
                                {% else %}
                                    <span class="not-available">Not Available</span>
                                {% endif %}
                            </td>
                            <td class="col-depth" {% if variant.Depth < 20 %}class="depth-warning"{% endif %}>
                                {{ variant.Depth }}
                            </td>
                            <td class="col-disease">{{ variant.Disease }}</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
            
            <div class="footer">
                <p>This report contains the top 300 prioritized variants. For full results, please refer to the complete Excel report.</p>
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
                {name: 'Chr', class: 'col-chr'},
                {name: 'Position', class: 'col-pos'},
                {name: 'Ref', class: 'col-ref'},
                {name: 'Alt', class: 'col-alt'},
                {name: 'Gene', class: 'col-gene'},
                {name: 'Location', class: 'col-location'},
                {name: 'dbSNP ID', class: 'col-dbsnp'},
                {name: 'COSMIC', class: 'col-cosmic'},
                {name: 'Transcript', class: 'col-transcript'},
                {name: 'CLNHGVS', class: 'col-clnhgvs'},
                {name: 'Filter', class: 'col-filter'},
                {name: 'Clinical Significance', class: 'col-significance'},
                {name: 'ACMG Classification', class: 'col-acmg-class'},
                {name: 'ACMG Met Criteria', class: 'col-acmg-criteria'},
                {name: 'Database', class: 'col-database'},
                {name: 'Sample AF', class: 'col-af'},
                {name: 'Population AF', class: 'col-gnomad-af'},
                {name: 'Depth', class: 'col-depth'},
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
                    
                    if (column === 'pos' || column === 'depth') {
                        const aVal = parseInt(a.querySelector(`.col-${column}`).textContent);
                        const bVal = parseInt(b.querySelector(`.col-${column}`).textContent);
                        return currentSort.ascending ? aVal - bVal : bVal - aVal;
                    }
                    
                    if (column === 'af' || column === 'gnomad') {
                        const aText = a.querySelector(`.col-${column}`).textContent;
                        const bText = b.querySelector(`.col-${column}`).textContent;
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
        elif 'likely_pathogenic' in significance:
            return 'likely-pathogenic'
        elif 'uncertain_significance' in significance:
            return 'uncertain'
        elif 'likely_benign' in significance:
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

    def clean_transcript_notation(self, transcript):
        """Clean transcript notation by removing special characters."""
        if pd.isna(transcript):
            return ""
        return str(transcript).replace('*', '')

    def get_database_links(self, gene, variant):
        """Generate database links for a given gene and variant."""
        chrom = str(variant.get('CHROM', '')).replace('chr', '')
        pos = str(variant.get('POS', ''))
        ref = str(variant.get('REF', ''))
        alt = str(variant.get('ALT', ''))
        variant_str = f"chr{chrom}-{pos}-{ref}-{alt}-hg19"
        hgvs_g = f"{variant.get('CHROM')}:g.{pos}{ref}>{alt}"
        
        # Determine ClinVar term priority: CLNHGVS > dbSNP ID > transcript
        clinvar_term = None
        if variant.get('CLNHGVS') and not pd.isna(variant.get('CLNHGVS')):
            clinvar_term = variant['CLNHGVS']
        elif variant.get('dbSNP_ID') and not pd.isna(variant.get('dbSNP_ID')) and variant['dbSNP_ID'] != 'Novel':
            clinvar_term = variant['dbSNP_ID']
        else:
            transcript = self.clean_transcript_notation(variant.get('Transcript', ''))
            if transcript:
                clinvar_term = transcript

        # Create database links including COSMIC
        links = {
            'ClinVar': f"https://www.ncbi.nlm.nih.gov/clinvar/?term={clinvar_term}" if clinvar_term else "#",
            'Franklin': f"https://franklin.genoox.com/clinical-db/variant/snp/{variant_str}",
            'Varsome': f"https://varsome.com/variant/hg19/{variant.get('dbSNP_ID') if variant.get('dbSNP_ID') and variant.get('dbSNP_ID') != 'Novel' else variant_str}?annotation-mode=germline",
            'OncoKB': f"https://www.oncokb.org/hgvsg/{hgvs_g}?refGenome=GRCh37"
        }
        
        # Add COSMIC link if COSMIC_ID exists
        if variant.get('COSMIC_ID') and not pd.isna(variant.get('COSMIC_ID')):
            cosmic_id = str(variant['COSMIC_ID']).replace('COSM', '')  # Remove 'COSM' prefix if present
            links['COSMIC'] = f"https://cancer.sanger.ac.uk/cosmic/mutation/overview?id={cosmic_id}"
        
        return links

    def get_filter_class(self, filter_value):
        """Determine CSS class based on FILTER value."""
        if pd.isna(filter_value):
            return 'filter-fail'
        return 'filter-pass' if filter_value == 'PASS' else 'filter-fail'

    def format_af(self, af_value):
        """Format allele frequency value."""
        try:
            if pd.isna(af_value):
                return "N/A"
            af = float(af_value)
            if af == 0:
                return "0.0000"
            return f"{af:.4f}"
        except (ValueError, TypeError):
            return "N/A"

    def count_variants_by_significance(self, df):
        """Count variants by their clinical significance."""
        pathogenic_count = 0
        benign_count = 0
        vus_count = 0
        conflicting_count = 0

        for significance in df['Clinical_Significance']:
            if pd.isna(significance):
                continue
            
            significance = str(significance).lower()
            if 'conflicting' in significance:
                conflicting_count += 1
            elif ('pathogenic' in significance and 'likely' not in significance) or \
                 'likely_pathogenic' in significance:
                pathogenic_count += 1
            elif 'uncertain_significance' in significance:
                vus_count += 1
            elif 'benign' in significance or 'likely_benign' in significance:
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
            f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Sample Allele Frequency">\n')
            f.write('##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="Population Allele Frequency (gnomAD)">\n')
            f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
            f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Symbol">\n')
            f.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical Significance">\n')
            f.write('##INFO=<ID=ACMG,Number=.,Type=String,Description="ACMG Classification">\n')
            f.write('##INFO=<ID=COSMIC,Number=1,Type=String,Description="COSMIC ID">\n')
            
            # Write VCF column headers
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            
            # Write variant data
            for _, row in df.iterrows():
                # Prepare INFO field
                info_fields = []
                if not pd.isna(row['AF']):
                    info_fields.append(f"AF={row['AF']}")
                if not pd.isna(row.get('gnomAD_AF')):
                    info_fields.append(f"gnomAD_AF={row['gnomAD_AF']}")
                if not pd.isna(row['Depth']):
                    info_fields.append(f"DP={int(row['Depth'])}")
                if not pd.isna(row['Gene']):
                    info_fields.append(f"GENE={row['Gene']}")
                if not pd.isna(row['Clinical_Significance']):
                    info_fields.append(f"CLNSIG={row['Clinical_Significance']}")
                if not pd.isna(row['ACMG_Classification']):
                    info_fields.append(f"ACMG={row['ACMG_Classification']}")
                if not pd.isna(row.get('COSMIC_ID')):
                    info_fields.append(f"COSMIC={row['COSMIC_ID']}")
                
                info = ';'.join(info_fields)
                
                # Write variant line
                variant_line = [
                    str(row['CHROM']),
                    str(row['POS']),
                    str(row['dbSNP_ID']) if not pd.isna(row['dbSNP_ID']) else '.',
                    str(row['REF']),
                    str(row['ALT']),
                    '.',  # QUAL
                    str(row['FILTER']) if not pd.isna(row['FILTER']) else '.',
                    info
                ]
                f.write('\t'.join(variant_line) + '\n')

    def generate_report(self, excel_file, output_file):
        """Generate HTML report and VCF from Excel file."""
        # Extract sample ID from filename
        sample_id = os.path.splitext(os.path.basename(excel_file))[0]
        
        # Read Excel file
        df = pd.read_excel(excel_file)
        
        # Select columns including COSMIC_ID
        columns = ['CHROM', 'POS', 'REF', 'ALT', 'Gene', 'Location', 'dbSNP_ID', 
                  'COSMIC_ID', 'Transcript', 'CLNHGVS', 'FILTER', 'Clinical_Significance', 
                  'ACMG_Classification', 'ACMG_Met_Criteria', 'Disease', 'AF', 'gnomAD_AF', 'Depth']
        df = df[columns]
        
        # Limit to top 300 variants
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
            
            # Add classes and formatting
            variant_dict['significance_class'] = self.get_significance_class(row['Clinical_Significance'])
            variant_dict['acmg_class'] = self.get_acmg_class(row['ACMG_Classification'])
            variant_dict['filter_class'] = self.get_filter_class(row['FILTER'])
            variant_dict['database_links'] = self.get_database_links(row['Gene'], variant_dict)
            
            # Format AF values
            variant_dict['AF_original'] = variant_dict['AF']
            variant_dict['AF_formatted'] = self.format_af(variant_dict['AF'])
            variant_dict['gnomAD_AF_formatted'] = self.format_af(variant_dict.get('gnomAD_AF'))
            
            # Convert Depth to int if possible
            try:
                variant_dict['Depth'] = int(variant_dict['Depth'])
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
            'pass_count': df['FILTER'].str.contains('PASS', case=True, na=False).sum(),
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
    if len(sys.argv) != 3:
        print("Usage: python report_generator.py input.xlsx output.html")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist")
        sys.exit(1)
    
    if not output_file.lower().endswith('.html'):
        output_file += '.html'
    
    try:
        generator = ClinicalReportGenerator()
        generator.generate_report(input_file, output_file)
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
