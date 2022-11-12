% %  Matlab package for analysis for fitting to bond lifetime vs force data %%
% %  -------------------------------------------------------------------------
% %  Georgia Tech license was used to use MATLAB
% %  ==============================================
% %  ----------------------------------------------
% %  Introduction of basic tutorial code
% %  ----------------------------------------------
% %  1. This code is simplified and summarized fitting.
% %    - There will be two boundaries when determining n* from fitting.
% %    - One is from directly applying Eq.S6 (for both)in the manuscript (mathmatically, minimum bound of n*)
% %    - The other is to consider one additional factor from molecular extension (maximum bound of n*) without any partial unfolding (i.e., n* = 0/ related to Eq.S7 for class 1, Eq. S13 for class2)
% %  2. For simplifying, this demo is implemented through scanning only discrete mean values of bond lifetime vs force data set
% %  3. One can calculate error by considering sem of bond lifetime of bond with Jacobian and residual matrix
% %  4. This script is encoded in "chronological order for analyzing"
% %  5. All sub-codes were written using Matlab-bulit in function or home-built codes
% %  6. One can select one of three Models

Any modification of code is available upon user's preference