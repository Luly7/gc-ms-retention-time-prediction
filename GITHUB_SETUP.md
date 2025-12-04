# GitHub Setup Guide

## Step 1: Initialize Git Repository

```bash
cd /Users/luly/Desktop/CAP10/Latest_Rt
git init
```

## Step 2: Configure Git (if not already done)

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

## Step 3: Add Files to Git

```bash
# Add all files (respecting .gitignore)
git add .

# Check what will be committed
git status
```

## Step 4: Create Initial Commit

```bash
git commit -m "Initial commit: GC-MS Retention Time Prediction System with 100 molecules and 20 column conditions"
```

## Step 5: Create GitHub Repository

### Option A: Using GitHub Website
1. Go to https://github.com
2. Click the "+" icon in top right corner
3. Select "New repository"
4. Fill in:
   - **Repository name**: `gc-ms-retention-time-prediction` (or your preferred name)
   - **Description**: "AI-powered GC-MS retention time prediction system using CrewAI with batch processing for 100 molecules across 20 column conditions"
   - **Visibility**: Choose Public or Private
   - **Do NOT** initialize with README, .gitignore, or license (we already have these)
5. Click "Create repository"

### Option B: Using GitHub CLI (if installed)
```bash
gh repo create gc-ms-retention-time-prediction --public --description "AI-powered GC-MS retention time prediction system"
```

## Step 6: Link Local Repository to GitHub

After creating the GitHub repo, you'll see commands like these (replace with your actual username):

```bash
git remote add origin https://github.com/YOUR_USERNAME/gc-ms-retention-time-prediction.git
git branch -M main
git push -u origin main
```

## Step 7: Verify Upload

```bash
# Check remote connection
git remote -v

# View commit history
git log --oneline
```

## Important: Protect Sensitive Information

### Before pushing, ensure your .env file is NOT committed:

```bash
# Check if .env exists in git
git ls-files | grep .env

# If it shows up, remove it:
git rm --cached .env
git commit -m "Remove .env from git"
```

### Create .env.example for others:

```bash
# Create a template
cat > .env.example << 'EOF'
# API Keys (replace with your actual keys)
OPENAI_API_KEY=your_openai_api_key_here
# ANTHROPIC_API_KEY=your_anthropic_key_here (optional)

# Optional settings
# CREWAI_TELEMETRY_OPT_OUT=true
EOF

git add .env.example
git commit -m "Add .env.example template"
git push
```

## Step 8: Update README for GitHub

The project already has comprehensive documentation:
- `README_BATCH_PROCESSING.md` - Technical documentation
- `QUICK_START.md` - User guide
- `EXPANSION_SUMMARY.md` - What was built

You may want to create a main README.md for GitHub visitors:

```bash
cat > README.md << 'EOF'
# GC-MS Retention Time Prediction System

AI-powered retention time prediction for Gas Chromatography-Mass Spectrometry (GC-MS) using CrewAI multi-agent system.

## Features

- 🧪 **100 Molecules Database** - Diverse chemical compounds across 11 classes
- 🔬 **20 Column Configurations** - Multiple GC column types and conditions
- 🤖 **AI-Powered Predictions** - Multi-agent system using CrewAI
- 📊 **Batch Processing** - Process multiple compounds simultaneously
- 📈 **Statistical Validation** - Comprehensive error analysis and reporting
- 📄 **Multiple Output Formats** - Markdown reports, CSV exports, visualizations

## Quick Start

\`\`\`bash
# Install dependencies
uv sync

# Set your API key
export OPENAI_API_KEY="your-key-here"

# Run the system
crewai run
\`\`\`

## Documentation

- **[Quick Start Guide](QUICK_START.md)** - Get started in 2 minutes
- **[Technical Documentation](README_BATCH_PROCESSING.md)** - Complete technical details
- **[Expansion Summary](EXPANSION_SUMMARY.md)** - System capabilities overview

## System Capabilities

- Process up to 100 compounds per batch
- Support for 20+ GC column configurations
- Generate 2,000+ predictions per run
- RDKit integration for molecular descriptors
- Statistical validation with multiple metrics
- CSV export for external analysis

## Requirements

- Python 3.10+
- OpenAI API key (or Anthropic API key)
- CrewAI
- Optional: RDKit for enhanced predictions

## Project Structure

\`\`\`
Latest_Rt/
├── knowledge/              # Data files
│   ├── molecules_database.csv (100 compounds)
│   └── column_conditions.csv (20 configs)
├── src/Latest_Rt/
│   ├── tools/             # Batch processing tools
│   ├── config/            # Agent and task configurations
│   ├── crew.py            # CrewAI setup
│   └── main.py            # Entry point
├── output/                # Generated reports
└── README.md              # This file
\`\`\`

## License

[Add your license here]

## Citation

If you use this system in your research, please cite appropriately.
EOF

git add README.md
git commit -m "Add main README for GitHub"
git push
```

## Future Updates

### To push changes later:

```bash
# After making changes
git add .
git commit -m "Describe your changes"
git push
```

### To pull changes from GitHub:

```bash
git pull origin main
```

### Create a new branch for features:

```bash
git checkout -b feature/new-feature
# Make changes
git add .
git commit -m "Add new feature"
git push -u origin feature/new-feature
# Then create Pull Request on GitHub
```

## Troubleshooting

### If push is rejected:
```bash
git pull origin main --rebase
git push
```

### If you need to change remote URL:
```bash
git remote set-url origin NEW_URL
```

### To view all branches:
```bash
git branch -a
```

## Additional Resources

- [GitHub Docs](https://docs.github.com)
- [Git Cheat Sheet](https://education.github.com/git-cheat-sheet-education.pdf)
- [CrewAI Documentation](https://docs.crewai.com)

