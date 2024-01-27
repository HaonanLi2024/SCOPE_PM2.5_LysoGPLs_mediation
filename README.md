# SCOPE_PM2.5_LysoGPLs_mediation

\<Lysoglycerophospholipids mediate the pro-atherosclerotic effects of fine particulate matter exposure>

- Background and Aim:

Recent studies highlight that lipids may mediate the pro-atherosclerotic effect of ambient fine particulate matter (PM2.5) exposure. Glycerophospholipid metabolism is the most widely recognized biological pathway disrupted by air pollution exposure. Lysoglycerophospholipids (LysoGPLs) belong to glycerophospholipids, and are known pro-atherosclerotic mediators.&#x20;

This study aims to investigate whether LysoGPLs mediate pro-atherosclerotic effect of PM2.5 exposure.

- Dataset:&#x20;

SCOPE was a panel study conducted in Beijing between August 2013 and  February 2015, with the primary aim of investigating the potential adverse effect of air pollution exposure on cardiometabolic biomarkers. The study recruited 120 participants aged 50 to 65 without a history of cardiovascular disease, residing in communities nearby Peking University campus. At enrollment, baseline demographic, and socioeconomic parameters of the participants were collected.

During the study period, ambient air pollutants were monitored, and each participant underwent up to seven repeated clinic visits at Peking University Hospital on the campus. Blood samples were collected.&#x20;

The datasets used in the present analysis are as following:

(a) PM2.5: Daily ambient PM2.5 concentration during the study period;

(b) LysoGPLs: Concentration of 18 LysoGPLs from 5 subclasses: including cPA 16:0/18:1, LPA 16:0/18:0/18:1/18:2/20:3/20:4/22:6, LPC(O) 16:0/18:0/18:1, LPG 16:0/ 16:1/18:0/18:1 and LysoPS 18:0/18:1;

(c) Pro-atherosclerotic clinical/subclinical biomarkers:&#x20;

Concentration of cytokines and chemokines, including interleukin-8 (IL-8), monocyte chemoattractant protein-1 (MCP-1), soluble CD40 ligand (sCD40L), interferon-γ (IFN-γ);

non-high-density lipoprotein cholesterol level;

platelet traits level, including mean platelet volume (MPV), platelet distribution width (PDW), platelet-large cell rate (P-LCR), and concentrations of platelet activation biomarkers, including sCD40L, thromboxane B2 (TxB2);

(d) Others: Meteorological variables, including temperature and relative humidity (RH); season, the day of week; baseline information of the participants, such as age, sex, body mass index (BMI), education level, income, marital status, smoking status, secondhand smoke exposure, and alcohol consumption frequency; ID of each participant.

- &#x20;Data Analysis Methods:

(a) Data preparation

(b) Data overview

(c) Statistical Analysis 1_linear mixed effect (LME) models

LME models were utilized to estimate the associations between short-term PM2.5 exposure and pro-atherosclerotic clinical/subclinical biomarkers or LysoGPLs.

(d) Statistical Analysis 2_mediation analyses

To further investigate whether LysoGPLs mediate the pro-atherosclerotic responses triggered by PM2.5, we conducted mediation analyses using a well-established two-stage regression approach based on the R package "mediation".

(e) Visualization with ggplot2

- Main Expected Results:

Table 1 Characteristics of the study participants at baseline.

Table 2 Concentrations of pro-atherosclerotic biomarkers and lysoglycerophospholipids in the study participants.

Table 3 Levels of fine particulate matter and meteorological parameters during the preceding 1 to 30 days before the clinical visits.

Figure 1 Percent changes in pro-atherosclerotic biomarkers and lysoglycerophospholipids per 10 μg/m3 increment of average PM2.5 1-30 days before clinical visits.

Figure 2 The proportions mediated by lysoglycerophospholipids on the associations between 14-d or 30-d average PM2.5 and pro-atherosclerotic biomarkers.

