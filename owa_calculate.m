clc
clear
format long g

%% Load data

[n,p] = uigetfile('*.txt','Select VCI data');
VCI = dlmread([p n]);

[n,p] = uigetfile('*.txt','Select TCI data');
TCI = dlmread([p n]);

[n,p] = uigetfile('*.txt','Select SMCI data');
SMCI = dlmread([p n]);

[n,p] = uigetfile('*.txt','Select PCI data');
PCI = dlmread([p n]);


%% Caclulate SE

% k : index 
% q : windows - kernel
% i : location (pixel)

SE_VCI = std(reshape(VCI,20*36,1))/ sqrt(20*36);
rk_VCI = 1-SE_VCI;

SE_TCI = std(reshape(TCI,20*36,1))/ sqrt(20*36) ;
rk_TCI = 1-SE_TCI;

SE_SMCI = std(reshape(SMCI,20*36,1))/ sqrt(20*36);
rk_SMCI = 1-SE_SMCI;

SE_PCI = std(reshape(PCI,20*36,1))/ sqrt(20*36);
rk_PCI = 1-SE_PCI;

%%
w_smci = 0.0423;
w_tci = 0.7179;
w_vci= 0.0101;
w_pci = 0.2296;

% alfa = 0.8

landa_SMCI = 0.494;
landa_PCI = 0.284;
landa_TCI = 0.157;
landa_VCI = 0.066;

K_size = input('enter Kernel size  : ');
St = (K_size/2)+0.5;
En = (K_size/2)-0.5;

for i = St:(size(VCI,1)-En)
    for j = St:(size(VCI,2)-En)
        
        % calculate V(alfa_ik_q)
        V_ik_q_VCI = (VCI(i,j)- min(min(VCI(i-En:i+En,j-En:j+En)))) / (max(max(VCI(i-En:i+En,j-En:j+En)))- min(min(VCI(i-En:i+En,j-En:j+En))));
        V_ik_q_TCI = (TCI(i,j)- min(min(TCI(i-En:i+En,j-En:j+En)))) / (max(max(TCI(i-En:i+En,j-En:j+En)))- min(min(TCI(i-En:i+En,j-En:j+En))));
        V_ik_q_SMCI = (SMCI(i,j)- min(min(SMCI(i-En:i+En,j-En:j+En)))) / (max(max(SMCI(i-En:i+En,j-En:j+En)))- min(min(SMCI(i-En:i+En,j-En:j+En))));
        V_ik_q_PCI = (PCI(i,j)- min(min(PCI(i-En:i+En,j-En:j+En)))) / (max(max(PCI(i-En:i+En,j-En:j+En)))- min(min(PCI(i-En:i+En,j-En:j+En))));
        
        % calculate W_kq
        rr_vci = VCI(i-En:i+En,j-En:j+En);
        r_kq_VCI =1-( std(reshape(rr_vci,size(rr_vci,1)*size(rr_vci,2),1)) / sqrt(size(rr_vci,1)*size(rr_vci,2)));
        W_st_kq_VCI = (w_vci * r_kq_VCI) / (rk_VCI);
        
        rr_tci = TCI(i-En:i+En,j-En:j+En);
        r_kq_TCI =1-( std(reshape(rr_tci,size(rr_tci,1)*size(rr_tci,2),1)) / sqrt(size(rr_tci,1)*size(rr_tci,2)));
        W_st_kq_TCI = (w_tci * r_kq_TCI) / (rk_TCI);
        
        rr_smci = TCI(i-En:i+En,j-En:j+En);
        r_kq_SMCI =1-( std(reshape(rr_smci,size(rr_smci,1)*size(rr_smci,2),1)) / sqrt(size(rr_smci,1)*size(rr_smci,2)));
        W_st_kq_SMCI = (w_tci * r_kq_SMCI) / (rk_SMCI);
        
        rr_pci = PCI(i-En:i+En,j-En:j+En);
        r_kq_PCI =1-( std(reshape(rr_pci,size(rr_pci,1)*size(rr_pci,2),1)) / sqrt(size(rr_pci,1)*size(rr_pci,2)));
        W_st_kq_PCI = (w_pci * r_kq_PCI) / (rk_PCI);
        
        
        W_kq_VCI = W_st_kq_VCI / (W_st_kq_VCI + W_st_kq_TCI + W_st_kq_SMCI + W_st_kq_PCI );
        
        W_kq_TCI = W_st_kq_TCI / (W_st_kq_VCI + W_st_kq_TCI + W_st_kq_SMCI + W_st_kq_PCI );
        
        W_kq_SMCI = W_st_kq_SMCI / (W_st_kq_VCI + W_st_kq_TCI + W_st_kq_SMCI + W_st_kq_PCI );
        
        W_kq_PCI = W_st_kq_PCI /(W_st_kq_VCI + W_st_kq_TCI + W_st_kq_SMCI + W_st_kq_PCI );
        
        
        landa_kol = [landa_SMCI landa_PCI landa_TCI landa_VCI];
        W_kq_kol = [W_kq_SMCI W_kq_PCI W_kq_TCI W_kq_VCI];
        V_ik_q_kol = [V_ik_q_SMCI V_ik_q_PCI V_ik_q_TCI V_ik_q_VCI];
        
        OWA(i,j) = sum(landa_kol .* W_kq_kol .* V_ik_q_kol)/sum(landa_kol .* W_kq_kol);
    end
end

%%

OWAA = OWA(St:end,St:end);
OWAAA = reshape(OWAA,size(OWAA,1)*size(OWAA,2),1);

[n,p] = uigetfile('*.txt','Select SPEI data');
speii = dlmread([p n]);

 spei =speii(St:size(speii,1)-En,St:size(speii,2)-En) ;
 SPEI = reshape(spei,size(spei,1)*size(spei,2),1);
r = corr2(SPEI,OWAAA);

