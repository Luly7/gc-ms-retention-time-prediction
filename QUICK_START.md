# Quick Start Guide - 100 Molecule Batch Processing

## What's New

Your CrewAI system has been expanded to handle **100 molecules** across **20 different GC column conditions**!

## Quick Run

```bash
cd /Users/luly/Desktop/CAP10/Latest_Rt
crewai run
```

## What Happens

The system will automatically:

1. **Extract Literature Data** (100 compounds)
   - Simulates mining retention time data from scientific papers
   - Includes DOI references, confidence scores
   - Covers diverse compound classes

2. **Generate Predictions** (2,000 predictions = 100 compounds × 20 columns)
   - Predicts retention times for each compound-column combination
   - Uses molecular descriptors (RDKit when available)
   - Accounts for column type, temperature program, film thickness

3. **Validate & Compare**
   - Compares predicted vs literature values
   - Calculates error statistics
   - Generates performance visualizations

4. **Create Reports**
   - Comprehensive validation report
   - Performance by compound class
   - CSV export for further analysis

## Data Files

### Input Data
- `knowledge/molecules_database.csv` - 100 compounds with SMILES, properties
- `knowledge/column_conditions.csv` - 20 column configurations

### Output Files (in `output/` directory)
- `latest_rt_literature_extraction.md` - Literature data summary
- `latest_rt_predictions.md` - Prediction results
- `latest_rt_validation_report.md` - Comprehensive validation
- `validation_results.csv` - Exportable data

## Molecule Database Highlights

**100 Compounds Including:**
- 🔷 Aromatic Hydrocarbons: Benzene, Toluene, Xylenes, Styrene, Naphthalene
- 🔶 Aliphatics: n-Hexane through n-Decane, Cycloalkanes
- 🟢 Alcohols: Methanol through Hexanol, Phenols, Cresols
- 🔴 Ketones: Acetone, Cyclohexanone, Acetophenone, Benzophenone
- 🟡 Acids & Esters: Acetic through Hexanoic acids, various esters
- 🟣 Amines: Methylamine, Aniline, Pyridines
- 🔵 Heterocyclics: Furan, Thiophene, Quinoline, Isoquinoline
- ⚫ Halogenated: Chlorobenzene, Bromobenzene, Chloroform
- ⚪ Others: Nitro compounds, Ethers, Polycyclic aromatics

## Column Configurations

**20 Different Setups:**
- DB-5, HP-5MS, Rtx-5, ZB-5 (non-polar)
- DB-1, SPB-1, Rtx-1 (100% dimethyl)
- DB-WAX (polar, polyethylene glycol)
- DB-35, DB-17, Rtx-35 (intermediate polarity)
- DB-624, Rtx-624 (cyanopropyl)
- Various temperature programs, film thicknesses

## New Tools Available

### Batch Processing
- `batch_extract_literature_rt` - Extract data for multiple compounds
- `batch_predict_retention_times` - Predict across all conditions
- `batch_compare_rt_values` - Statistical validation

### Reporting
- `generate_summary_report` - Markdown/JSON reports
- `create_performance_visualization` - ASCII charts
- `export_results_to_csv` - Data export
- `generate_compound_class_report` - Class-specific analysis

## Performance Metrics

Expected accuracy distribution:
- ✅ Excellent (<5% error): 40-60% of predictions
- 👍 Good (5-10% error): 25-35% of predictions
- 🆗 Moderate (10-20% error): 10-20% of predictions
- ⚠️ Poor (>20% error): <10% of predictions

## Customization

### Process Fewer Compounds
Edit `src/Latest_Rt/config/tasks.yaml` and change `num_compounds` parameter:
```yaml
Use the batch_extract_literature_rt tool with num_compounds=50
```

### Select Specific Columns
Change `column_conditions` parameter:
```yaml
column_conditions="1,2,3,4"  # Only first 4 columns
```

## Troubleshooting

### "Module not found" errors
Make sure you're in the correct directory:
```bash
cd /Users/luly/Desktop/CAP10/Latest_Rt
```

### Slow processing
- Normal: 2-5 minutes for full run
- Depends on API response times
- Can reduce compound count if needed

### Want better predictions?
Install RDKit for real molecular descriptor calculation:
```bash
uv add rdkit
```

## Next Steps

1. **Review Output Files** - Check `output/` directory
2. **Analyze CSV Data** - Open in Excel or Python pandas
3. **Customize Parameters** - Edit YAML configs for your needs
4. **Add More Molecules** - Extend `molecules_database.csv`
5. **Add Column Types** - Extend `column_conditions.csv`

## Example Output Structure

```
output/
├── latest_rt_literature_extraction.md
│   └── Summary of 100 compounds with literature RT data
├── latest_rt_predictions.md
│   └── 2,000 predictions with analysis
├── latest_rt_validation_report.md
│   └── Comprehensive validation with statistics
└── validation_results.csv
    └── All results in spreadsheet format
```

## Tips for Best Results

1. **Let it run completely** - Don't interrupt the process
2. **Check all output files** - Each has different insights
3. **Use CSV for analysis** - Easy to import into other tools
4. **Review by compound class** - Some classes predict better than others
5. **Compare columns** - See which columns work best for your compounds

## Support

For detailed information, see:
- `README_BATCH_PROCESSING.md` - Complete technical documentation
- `knowledge/molecules_database.csv` - Full compound list
- `knowledge/column_conditions.csv` - All column specs

---

**Ready to analyze 100 molecules? Just run `crewai run`!** 🚀

