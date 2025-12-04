# GC-MS Retention Time Prediction - Batch Processing System

## Overview

This system has been expanded to handle **100 molecules** across **20 different column conditions**, providing comprehensive retention time predictions and validation for analytical chemistry workflows.

## System Capabilities

### 1. Molecule Database (100 Compounds)
Located: `knowledge/molecules_database.csv`

**Compound Classes Included:**
- Aromatic Hydrocarbons (Benzene, Toluene, Xylenes, etc.)
- Polycyclic Aromatics (Naphthalene, Anthracene, Phenanthrene)
- Aliphatic Hydrocarbons (n-Hexane through n-Decane)
- Cycloalkanes (Cyclohexane, Methylcyclohexane, etc.)
- Alcohols (Methanol through Hexanol, including branched)
- Phenols (Phenol, Cresols)
- Ketones (Acetone, Cyclohexanone, Aromatic ketones)
- Carboxylic Acids (Acetic through Hexanoic, Benzoic)
- Esters (Methyl/Ethyl acetates, propionates, benzoates)
- Aldehydes (Acetaldehyde, Benzaldehyde, Cinnamaldehyde)
- Amines (Methylamine, Aniline, Dimethylaniline)
- Heterocyclics (Pyridine, Quinoline, Furan, Thiophene)
- Halogenated Compounds (Chlorobenzene, Bromobenzene, etc.)
- Nitro Compounds (Nitrobenzene, Nitrotoluenes)
- Ethers (Anisole, Diethyl ether, THF, Dioxane)

**Data Fields:**
- Compound ID, Name, SMILES notation
- Molecular weight, CAS number
- Compound class classification

### 2. Column Conditions Database (20 Configurations)
Located: `knowledge/column_conditions.csv`

**Column Types:**
- **DB-5** (5% Phenyl, 95% Dimethylpolysiloxane) - 3 configurations
- **HP-5MS** (5% Phenyl, MS-optimized) - 2 configurations
- **Rtx-5** (5% Diphenyl) - 1 configuration
- **ZB-5** (5% Phenyl) - 1 configuration
- **DB-1** (100% Dimethylpolysiloxane) - 2 configurations
- **DB-WAX** (Polyethylene Glycol) - 2 configurations
- **SPB-1** (100% Dimethylpolysiloxane) - 1 configuration
- **Rtx-1** (100% Dimethylpolysiloxane) - 1 configuration
- **DB-35** (35% Phenyl) - 1 configuration
- **DB-17** (50% Phenyl) - 1 configuration
- **Rtx-35** (35% Diphenyl) - 1 configuration
- **DB-624** (6% Cyanopropylphenyl) - 1 configuration
- **Rtx-624** (6% Cyanopropylphenyl) - 1 configuration
- **DB-VRX** (20% Phenyl) - 1 configuration
- **BPX5** (5% Phenyl Polysilphenylene) - 1 configuration

**Variable Parameters:**
- Column length: 20-30m
- Internal diameter: 0.18-0.25mm
- Film thickness: 0.25-1.40μm
- Temperature programs: Various ramp rates (8-15°C/min)
- Carrier gas: Helium or Nitrogen
- Flow rates: 1.0-1.5 mL/min

### 3. Batch Processing Tools

#### `batch_extract_literature_rt`
Extracts retention time data from literature for multiple compounds.
- **Input**: Number of compounds (default 10, max 100)
- **Output**: JSON with literature RT data, DOI references, confidence scores
- **Features**: Simulates literature mining with realistic RT values

#### `batch_predict_retention_times`
Predicts retention times for multiple compounds across column conditions.
- **Input**: Number of compounds, column conditions ("all" or specific IDs)
- **Output**: Comprehensive predictions for all compound-column combinations
- **Features**: 
  - Uses RDKit for molecular descriptor calculation (if available)
  - Column-specific factors (stationary phase effects)
  - Temperature program adjustments
  - Film thickness corrections
  - Generates up to 2,000 predictions (100 compounds × 20 columns)

#### `batch_compare_rt_values`
Compares experimental vs predicted RT with statistical analysis.
- **Input**: Literature data JSON, Prediction data JSON
- **Output**: Comprehensive validation statistics
- **Metrics**:
  - Mean/std deviation of relative errors
  - Agreement distribution (Excellent/Good/Moderate/Poor)
  - Per-compound comparison details
  - Performance by compound class

#### `generate_summary_report`
Creates formatted summary reports from comparison data.
- **Output formats**: Markdown or JSON
- **Includes**:
  - Executive summary with key metrics
  - Agreement distribution tables
  - Top 10 best/worst predictions
  - Recommendations for improvement

### 4. Visualization Tools

#### `create_performance_visualization`
Generates ASCII-based charts and graphs.
- Agreement distribution bar chart
- Error distribution histogram
- Performance by compound class table

#### `export_results_to_csv`
Exports validation results to CSV for external analysis.
- Compatible with Excel, R, Python pandas
- All compound details and error metrics included

#### `generate_compound_class_report`
Detailed analysis grouped by compound class.
- Average error per class
- Standard deviation and range
- Agreement distribution per class
- List of compounds in each class

## Usage

### Running the Batch Analysis

```bash
cd /Users/luly/Desktop/CAP10/Latest_Rt
crewai run
```

The system will automatically:
1. Extract literature RT data for 100 compounds
2. Generate predictions across 20 column conditions (2,000 predictions)
3. Compare and validate all predictions
4. Generate comprehensive reports with visualizations

### Output Files

All outputs are saved to `output/` directory:
- `latest_rt_literature_extraction.md` - Literature mining results
- `latest_rt_predictions.md` - Prediction results and analysis
- `latest_rt_validation_report.md` - Comprehensive validation report
- `validation_results.csv` - Exportable data for further analysis

## Performance Expectations

### Prediction Accuracy Targets
- **Excellent** (<5% error): Target for 40-60% of predictions
- **Good** (5-10% error): Target for 25-35% of predictions
- **Moderate** (10-20% error): Expected for 10-20% of predictions
- **Poor** (>20% error): Should be <10% of predictions

### Processing Capacity
- **100 compounds** × **20 column conditions** = **2,000 predictions**
- Processing time: ~2-5 minutes (depending on LLM response time)
- Memory usage: Minimal (<100MB)

## Technical Details

### Molecular Descriptor Calculation
When RDKit is available, the system calculates:
- Molecular weight, LogP, TPSA
- H-bond donors/acceptors
- Rotatable bonds
- Aromatic/aliphatic/saturated rings
- Molar refractivity
- Balaban J index
- Labute ASA
- PEOE VSA descriptors

### Column Factor Adjustments
- **DB-5/HP-5MS**: Baseline (factor 1.0-1.05)
- **DB-1**: Slightly faster (factor 0.95)
- **DB-WAX**: Significantly slower for polar compounds (factor 1.25)
- **DB-35/DB-17**: Intermediate polarity (factors 1.08-1.12)

### Temperature Program Effects
- Fast ramps (15°C/min): Shorter retention times (factor 0.95)
- Slow ramps (8°C/min): Longer retention times (factor 1.05)
- Standard (10°C/min): Baseline (factor 1.0)

## Future Enhancements

1. **Real Literature Integration**: Connect to PubChem, ChemSpider, NIST databases
2. **Machine Learning Models**: Train on real experimental data
3. **Interactive Visualizations**: Add plotly/matplotlib charts
4. **API Endpoints**: RESTful API for external integrations
5. **Database Backend**: PostgreSQL/MongoDB for large-scale data
6. **Web Interface**: Flask/FastAPI dashboard for easy access

## Dependencies

```bash
# Core dependencies (already installed)
crewai[anthropic]
pydantic
python-dotenv

# Optional for enhanced functionality
rdkit  # For molecular descriptor calculation
pandas  # For advanced data analysis
matplotlib  # For plotting
```

## Troubleshooting

### RDKit Not Available
If RDKit is not installed, the system falls back to simplified molecular weight-based predictions. To install:
```bash
uv add rdkit
```

### Memory Issues
For very large batch sizes (>100 compounds), consider processing in chunks by modifying the `num_compounds` parameter.

### API Rate Limits
If using Anthropic API, be aware of rate limits. The system processes sequentially to avoid issues.

## Citation

If you use this system in your research, please cite:
```
GC-MS Retention Time Prediction System
CrewAI-based Multi-Agent System for Analytical Chemistry
2024
```

## License

This project is part of the CAP10 research initiative.

## Contact

For questions or issues, please refer to the main project documentation.

