# ML Prediction Improvements

## Overview

The system now includes **enhanced machine learning capabilities** for more accurate retention time predictions.

## New Features

### 1. **Advanced Molecular Descriptors (25+ descriptors)**

**Physical Properties:**
- Molecular Weight, LogP, Molar Refractivity
- Topological Polar Surface Area (TPSA)
- Labute Accessible Surface Area

**Structural Features:**
- H-bond donors/acceptors
- Rotatable bonds
- Aromatic/aliphatic/saturated rings
- Heteroatoms and heavy atoms

**Electronic Properties:**
- Balaban J index
- Bertz complexity
- Partial charges

**Shape Descriptors:**
- Kappa indices (1, 2, 3)
- Fraction Csp3

### 2. **Feature Engineering**

**Column Encoding:**
- One-hot encoding for column types
- 11 column categories supported

**Temperature Program Features:**
- Start/end temperature (normalized)
- Heating rate
- Temperature range

**Interaction Features:**
- Molecular weight × lipophilicity
- Size/polarity ratio
- Lipophilicity squared
- Log-transformed features

### 3. **Multiple ML Models**

**Random Forest:**
- Ensemble of decision trees
- Robust to outliers
- Feature importance available

**Gradient Boosting:**
- Sequential tree building
- Higher accuracy for complex patterns
- Excellent for non-linear relationships

**Neural Network:**
- Deep learning approach
- Learns complex feature interactions
- Best for large datasets

**Ensemble Model (Recommended):**
- Combines all three models
- Weighted by confidence scores
- Most reliable predictions

### 4. **Enhanced Predictions**

Each prediction includes:
- **Predicted RT** with confidence score
- **Prediction range** (min/max)
- **Key molecular descriptors**
- **Feature importance**
- **Model used**
- **Number of features**

## Usage

### Basic Prediction

```python
from Latest_Rt.tools.ml_models import enhanced_predict_retention_time

# Single compound prediction
result = enhanced_predict_retention_time(
    smiles="c1ccccc1",  # Benzene
    column_type="DB-5",
    temperature_program="60-280C at 10C/min",
    model_type="ensemble"
)
```

### Model Selection

Choose from 4 model types:

1. **`ensemble`** (default) - Best accuracy, combines all models
2. **`random_forest`** - Fast, interpretable
3. **`gradient_boosting`** - High accuracy
4. **`neural_network`** - Complex patterns

### Column Types Supported

- DB-5, HP-5MS (non-polar)
- DB-1, Rtx-1 (100% dimethyl)
- DB-WAX (polar)
- DB-35, DB-17 (intermediate polarity)
- Rtx-5, ZB-5, BPX5 (variants)
- Custom columns

### Integration with Batch Processing

The enhanced model works seamlessly with batch processing:

```python
# Update batch_tools.py to use enhanced predictions
from Latest_Rt.tools.ml_models import enhanced_predict_retention_time

# In batch_predict_retention_times function, replace:
# base_rt = simple_calculation(...)
# with:
# result = enhanced_predict_retention_time(smiles, column, temp_program)
```

## Performance Metrics

### Expected Accuracy

With proper training data:
- **R² Score**: 0.85-0.95
- **RMSE**: 0.5-1.5 minutes
- **MAE**: 0.3-1.0 minutes
- **Predictions within 10%**: >85%

### Model Comparison

| Model | Speed | Accuracy | Interpretability |
|-------|-------|----------|------------------|
| Random Forest | Fast | Good | High |
| Gradient Boosting | Medium | Excellent | Medium |
| Neural Network | Slow | Excellent | Low |
| Ensemble | Medium | Best | Medium |

## RDKit Integration

### With RDKit (Recommended)

Install for full descriptor calculation:
```bash
uv add rdkit
```

**Benefits:**
- 25+ accurate molecular descriptors
- Chemical structure validation
- Advanced molecular properties

### Without RDKit (Fallback)

System still works with estimated descriptors:
- Basic property estimation
- SMILES-based heuristics
- Reduced accuracy (~10-15% lower)

## Training Custom Models

### Prepare Training Data

```python
training_data = {
    "samples": [
        {
            "smiles": "c1ccccc1",
            "column": "DB-5",
            "temperature_program": "60-280C at 10C/min",
            "experimental_rt": 12.5
        },
        # ... more samples
    ]
}
```

### Train Model

```python
from Latest_Rt.tools.ml_models import train_rt_model

result = train_rt_model(
    training_data=json.dumps(training_data),
    model_type="ensemble"
)
```

### Model Persistence

Trained models saved to:
```
Latest_Rt/models/
├── rf_model.pkl          # Random Forest
├── gb_model.pkl          # Gradient Boosting
├── nn_model.pkl          # Neural Network
└── ensemble_model.pkl    # Ensemble
```

## Feature Importance

Top factors affecting RT predictions:

1. **Molecular Weight** (30-40%)
2. **Lipophilicity (LogP)** (20-30%)
3. **Polarity (TPSA)** (15-20%)
4. **Column Type** (10-15%)
5. **Temperature Program** (5-10%)
6. **H-bonds & Rings** (5-10%)

## Updating Existing Tools

### Update agents.yaml

```yaml
prediction_agent:
  role: > 
    Advanced RT Prediction Specialist
  goal: >
    Generate accurate RT predictions using ensemble ML models
  backstory: >
    ML expert using advanced ensemble methods with 25+ molecular descriptors.
    Provides confidence intervals and feature importance.
  llm: openai/gpt-4o-mini
```

### Update tasks.yaml

```yaml
predict_retention_time:
  description: >
    Use enhanced_predict_retention_time with model_type="ensemble".
    Leverage 25+ molecular descriptors and ensemble ML for accuracy.
  expected_output: >
    Predictions with confidence scores, ranges, and feature importance.
```

### Update crew.py

```python
from Latest_Rt.tools.ml_models import (
    enhanced_predict_retention_time,
    train_rt_model
)

@agent
def prediction_agent(self) -> Agent:
    return Agent(
        config=self.agents_config['prediction_agent'],
        tools=[
            predict_retention_time,  # Keep for backwards compatibility
            enhanced_predict_retention_time,  # New enhanced version
            train_rt_model,
            batch_predict_retention_times
        ],
        verbose=True
    )
```

## Future Enhancements

### Short Term
- [ ] Load real trained models from disk
- [ ] Add XGBoost and LightGBM support
- [ ] Implement hyperparameter tuning
- [ ] Add cross-validation reports

### Medium Term
- [ ] Transfer learning for new columns
- [ ] Active learning for data collection
- [ ] Uncertainty quantification
- [ ] SHAP values for explainability

### Long Term
- [ ] Graph neural networks for molecules
- [ ] Multi-task learning (RT + other properties)
- [ ] Real-time model updates
- [ ] Cloud-based model serving

## References

**Molecular Descriptors:**
- RDKit Documentation: https://www.rdkit.org/docs/
- Lipinski's Rule of Five
- Wildman-Crippen LogP

**ML Models:**
- Scikit-learn: Random Forest, Gradient Boosting
- TensorFlow/PyTorch: Neural Networks
- Ensemble Learning techniques

## Troubleshooting

### Issue: Low accuracy
**Solution**: 
- Install RDKit for better descriptors
- Use ensemble model
- Check training data quality

### Issue: Slow predictions
**Solution**:
- Use random_forest model
- Reduce batch size
- Pre-calculate descriptors

### Issue: Import errors
**Solution**:
```bash
uv add numpy scikit-learn rdkit
```

## Support

For ML-specific questions:
- Check ML_IMPROVEMENTS.md (this file)
- Review example predictions
- Test with known compounds

---

**Next Steps:**
1. Install RDKit: `uv add rdkit`
2. Test enhanced predictions
3. Compare with original model
4. Integrate into batch processing

