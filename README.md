# Longitudinal Analysis and Prediction of Functional Decline in the ALFA+ Cohort

This repository contains the complete analytical pipeline for my 2026  Master's Thesis. The project investigates the longitudinal trajectories of biomarkers associated with functional decline and develops predictive models for clinically meaningful decline using data from the ALFA+ cohort.

From an analytical and reproducibility perspective, all work is documented and structured into three dedicated R scripts, each corresponding to the main milestones of the project.

---

## 📋 Table of Contents
- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)

---

## 🔬 Project Overview

The objective of this study is twofold:
1. **Longitudinal Assessment:** To model how specific biomarkers influence the trajectory of functional decline over time using Linear Mixed-Effects Models (LMMs).
2. **Predictive Modeling:** To classify and predict whether a patient will experience a clinically meaningful decline using logistic regression and evaluate performance via Receiver Operating Characteristic (ROC) - Area Under the Curve (AUC) methodologies.

Functional decline is operationalized throughout using the **Amsterdam Instrumental Activities of Daily Living Questionnaire (A-IADL-Q)**.

---

## 📁 Repository Structure

```text
├── README.md
├── ALFA+_dataset_create_program.R
├── ALFA+_A-IADL-Q_longitudinal_analysis.R
└── ALFA+_clinically_meaningful_decline_prediction.R
