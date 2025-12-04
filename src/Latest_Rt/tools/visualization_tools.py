from crewai.tools import tool
import json
from pathlib import Path

@tool("Create Performance Visualization")
def create_performance_visualization(comparison_data: str) -> str:
    """
    Create ASCII-based visualizations of model performance.
    Generates charts and graphs for validation results.
    
    Args:
        comparison_data: JSON string from batch_compare_rt_values
    
    Returns:
        String with ASCII visualization charts
    """
    try:
        data = json.loads(comparison_data)
        
        if data.get('status') == 'error':
            return "Error: Unable to create visualization from error data"
        
        stats = data['summary_statistics']
        dist = stats['agreement_distribution']
        
        # Create bar chart for agreement distribution
        chart = """
## Performance Visualization

### Agreement Distribution Bar Chart

"""
        
        max_count = max(dist['excellent_count'], dist['good_count'], 
                       dist['moderate_count'], dist['poor_count'])
        scale = 50 / max_count if max_count > 0 else 1
        
        categories = [
            ("Excellent (<5%)", dist['excellent_count'], dist['excellent_percent']),
            ("Good (5-10%)", dist['good_count'], dist['good_percent']),
            ("Moderate (10-20%)", dist['moderate_count'], dist['moderate_percent']),
            ("Poor (>20%)", dist['poor_count'], dist['poor_percent'])
        ]
        
        for label, count, percent in categories:
            bar_length = int(count * scale)
            bar = "█" * bar_length
            chart += f"{label:20} | {bar} {count} ({percent}%)\n"
        
        # Error distribution histogram
        chart += """
### Relative Error Distribution

"""
        
        # Group errors into bins
        comparisons = data['detailed_comparisons']
        bins = [0, 5, 10, 15, 20, 25, 30, 100]
        bin_counts = [0] * (len(bins) - 1)
        
        for comp in comparisons:
            error = comp['relative_error_percent']
            for i in range(len(bins) - 1):
                if bins[i] <= error < bins[i + 1]:
                    bin_counts[i] += 1
                    break
        
        max_bin = max(bin_counts) if bin_counts else 1
        scale = 40 / max_bin if max_bin > 0 else 1
        
        for i in range(len(bins) - 1):
            if bins[i + 1] == 100:
                label = f"{bins[i]}%+"
            else:
                label = f"{bins[i]}-{bins[i+1]}%"
            bar_length = int(bin_counts[i] * scale)
            bar = "▓" * bar_length
            chart += f"{label:10} | {bar} {bin_counts[i]}\n"
        
        # Compound class performance
        chart += """
### Performance by Compound Class

"""
        
        class_stats = {}
        for comp in comparisons:
            comp_class = comp['compound_class']
            if comp_class not in class_stats:
                class_stats[comp_class] = {'count': 0, 'total_error': 0, 'excellent': 0}
            
            class_stats[comp_class]['count'] += 1
            class_stats[comp_class]['total_error'] += comp['relative_error_percent']
            if comp['agreement_level'] == 'Excellent':
                class_stats[comp_class]['excellent'] += 1
        
        # Sort by average error
        sorted_classes = sorted(class_stats.items(), 
                               key=lambda x: x[1]['total_error'] / x[1]['count'])
        
        chart += "| Compound Class | Avg Error | Count | Excellent % |\n"
        chart += "|----------------|-----------|-------|-------------|\n"
        
        for comp_class, stats_dict in sorted_classes[:10]:  # Top 10 classes
            avg_error = stats_dict['total_error'] / stats_dict['count']
            excellent_pct = (stats_dict['excellent'] / stats_dict['count'] * 100)
            chart += f"| {comp_class:25} | {avg_error:5.2f}% | {stats_dict['count']:5} | {excellent_pct:6.1f}% |\n"
        
        return chart
        
    except Exception as e:
        return f"Error creating visualization: {str(e)}"


@tool("Export Results to CSV")
def export_results_to_csv(comparison_data: str, output_filename: str = "validation_results.csv") -> str:
    """
    Export validation results to CSV format for further analysis.
    
    Args:
        comparison_data: JSON string from batch_compare_rt_values
        output_filename: Name of the output CSV file
    
    Returns:
        Status message with file location
    """
    try:
        import csv
        
        data = json.loads(comparison_data)
        
        if data.get('status') == 'error':
            return "Error: Unable to export error data"
        
        # Get output directory
        output_dir = Path(__file__).parent.parent.parent.parent / "output"
        output_dir.mkdir(exist_ok=True)
        
        output_path = output_dir / output_filename
        
        # Write CSV
        with open(output_path, 'w', newline='') as csvfile:
            fieldnames = [
                'compound_id', 'compound_name', 'compound_class',
                'experimental_rt', 'predicted_rt', 'absolute_error',
                'relative_error_percent', 'agreement_level', 'column',
                'literature_source'
            ]
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for comp in data['detailed_comparisons']:
                writer.writerow(comp)
        
        return f"Successfully exported {len(data['detailed_comparisons'])} results to {output_path}"
        
    except Exception as e:
        return f"Error exporting to CSV: {str(e)}"


@tool("Generate Compound Class Report")
def generate_compound_class_report(comparison_data: str) -> str:
    """
    Generate detailed performance analysis grouped by compound class.
    
    Args:
        comparison_data: JSON string from batch_compare_rt_values
    
    Returns:
        Markdown report analyzing performance by compound class
    """
    try:
        data = json.loads(comparison_data)
        
        if data.get('status') == 'error':
            return "Error: Unable to generate report from error data"
        
        # Group by compound class
        class_data = {}
        for comp in data['detailed_comparisons']:
            comp_class = comp['compound_class']
            if comp_class not in class_data:
                class_data[comp_class] = {
                    'compounds': [],
                    'errors': [],
                    'excellent': 0,
                    'good': 0,
                    'moderate': 0,
                    'poor': 0
                }
            
            class_data[comp_class]['compounds'].append(comp['compound_name'])
            class_data[comp_class]['errors'].append(comp['relative_error_percent'])
            
            agreement = comp['agreement_level']
            if agreement == 'Excellent':
                class_data[comp_class]['excellent'] += 1
            elif agreement == 'Good':
                class_data[comp_class]['good'] += 1
            elif agreement == 'Moderate':
                class_data[comp_class]['moderate'] += 1
            else:
                class_data[comp_class]['poor'] += 1
        
        # Generate report
        report = "# Performance Analysis by Compound Class\n\n"
        
        # Sort by average error
        sorted_classes = sorted(class_data.items(),
                               key=lambda x: sum(x[1]['errors']) / len(x[1]['errors']))
        
        for comp_class, stats in sorted_classes:
            errors = stats['errors']
            n_compounds = len(errors)
            avg_error = sum(errors) / n_compounds
            min_error = min(errors)
            max_error = max(errors)
            
            # Calculate standard deviation
            variance = sum((e - avg_error) ** 2 for e in errors) / n_compounds
            std_dev = variance ** 0.5
            
            report += f"## {comp_class}\n\n"
            report += f"**Number of Compounds**: {n_compounds}\n\n"
            report += f"**Performance Metrics**:\n"
            report += f"- Average Relative Error: {avg_error:.2f}%\n"
            report += f"- Standard Deviation: {std_dev:.2f}%\n"
            report += f"- Best Prediction: {min_error:.2f}%\n"
            report += f"- Worst Prediction: {max_error:.2f}%\n\n"
            
            report += f"**Agreement Distribution**:\n"
            report += f"- Excellent (<5%): {stats['excellent']} ({stats['excellent']/n_compounds*100:.1f}%)\n"
            report += f"- Good (5-10%): {stats['good']} ({stats['good']/n_compounds*100:.1f}%)\n"
            report += f"- Moderate (10-20%): {stats['moderate']} ({stats['moderate']/n_compounds*100:.1f}%)\n"
            report += f"- Poor (>20%): {stats['poor']} ({stats['poor']/n_compounds*100:.1f}%)\n\n"
            
            # List compounds
            report += f"**Compounds in this class**: {', '.join(stats['compounds'][:5])}"
            if n_compounds > 5:
                report += f" and {n_compounds - 5} more"
            report += "\n\n"
            
            report += "---\n\n"
        
        return report
        
    except Exception as e:
        return f"Error generating compound class report: {str(e)}"

