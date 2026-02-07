#!/bin/bash

# 1. Install Python dependencies
echo "Installing Python libraries..."
pip install -r requirements.txt

# 2. Install R dependencies (NEW)
echo "Installing R libraries (this may take a moment)..."
Rscript -e "install.packages(c('jsonlite', 'ordinal', 'insight'), repos='https://cloud.r-project.org')"

# 3. Authenticate with Google (User Login)
echo "Starting Google Authentication..."
gcloud auth login

# 4. Set the Quota Project (Vital for gspread 403 errors)
echo "Setting Quota Project to 'blues-medical'..."
gcloud auth application-default set-quota-project blues-medical

# 5. Create Application Default Credentials with Drive/Sheets scopes
echo "Authorizing Drive and Sheets access..."
gcloud auth application-default login --scopes=https://www.googleapis.com/auth/drive,https://www.googleapis.com/auth/spreadsheets,https://www.googleapis.com/auth/cloud-platform

echo "---------------------------------------"
echo "Setup Complete! You can now run your scripts."