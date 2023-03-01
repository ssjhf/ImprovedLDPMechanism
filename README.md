# ImprovedLDPMechanism
1.This repository is the supplemental material for the paper---"An Improved Christofides Mechanism for Local Differential Privacy Framework" which was submitted to VLDB2023.

2.Dataset
The health_insurance.dat is the HCOVANY dataset, which indicates whether persons had any health insurance coverage during the interview, as measured by employer-provided insurance(HINSEMP), privately purchased insurance (HINSPUR), Medicare (HINSCARE), Medicaid or other governmental insurance (HINSCAID), TRICARE or other military care (HINSTRI), or Veterans Administration-provided insurance (HINSVA). The Census Bureau does not consider the respondents with health insurance coverage if their only coverage is from Indian Health Services (HINSIHS), as IHS policies are not always comprehensive. Codes of 2 indicate that a person is covered (either directly or through another household member's policy) by the given type of insurance; codes of 1 indicate that a person is not covered. We select the subset of HCOVANY Dataset in our experiment (the health insurance coverage of 1% of total population in USA in 2021, the sampling size N=3252599).

3.Source Code
Fig1.m and Fig5.m are the MATLAB programs used to plot figure 1 and 5, respectively. Fig234.m is the MATLAB program used to plot figure 2 to 4. Table1.m and Table2.m are the MATLAB programs used to get Table 1 and 2, respectively.
