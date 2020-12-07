% T = table(H1_C1.Position,H2_C1.Position,Hv_C1.Position,H1_C2.Position,H2_C2.Position,Hv_C2.Position,H1_C3.Position,H2_C3.Position,Hv_C3.Position);
H1_C1 = H1_C1.Position'
T1 = table(H1_C1.Position);
filename = 'amplitude.xlsx';
writetable(T1,filename);