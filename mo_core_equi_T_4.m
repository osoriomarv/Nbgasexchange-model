function [dIW,A]= mo_core_equi_T_4(SiO2, MgO, FeO, Al2O3, CaO, NiO, Fe, Ni, O, Si, Tadi, Padi, frac)
% Take D value as .01, explore D value, show the change at low Dhe - How
% solution varies with core size and the influence of IW. 
% Use enstatite as nominal compositon

melt_mols_i = [SiO2, MgO, FeO, Al2O3*2, CaO, NiO];
metal_mols_i_scaled=[Fe, Ni, O, Si]*frac;
P=Padi;
T=Tadi;

FeO_guess_v=logspace(-1.5,-.01,1e6);

% 
% melt_mols_i = [4.762, 4.3612, 0, .45031, 0.21209, 0];
% metal_mols_i_scaled = [5.1929, .29816, 1e-10, 1.1901];
% 
% %Pass into mass of core and bse
% P = 60;
% T = 1833;


T_ref=1873;
o=0;

%partitioning expressios
KD_Si_p=[1.3,-13500,0];
KD_Ni_p=[0.46,2700,-61];
KD_O_p=[0.6,-3800,22];

E_O_O_ref=0;
E_O_Si_ref=0;
E_Si_Si_ref=0;


melt_mol_frac_i=melt_mols_i./sum(melt_mols_i);
metal_mol_frac_i=metal_mols_i_scaled./sum(metal_mols_i_scaled);


FeO_i=melt_mols_i(3);
NiO_i=melt_mols_i(6);
SiO2_i=melt_mols_i(1);
MgO_i=melt_mols_i(2);
AlO1_5_i=melt_mols_i(4);
CaO_i=melt_mols_i(5);

Fe_i=metal_mols_i_scaled(1);
Ni_i=metal_mols_i_scaled(2);
O_i=metal_mols_i_scaled(3);
Si_i=metal_mols_i_scaled(4);

Fe_mol_frac_pick(1)=metal_mol_frac_i(1);
Ni_mol_frac_pick(1)=metal_mol_frac_i(2);
O_mol_frac_pick(1)=metal_mol_frac_i(3);
Si_mol_frac_pick(1)=metal_mol_frac_i(4);

    E_O_O=E_O_O_ref.*T_ref./T;
    E_O_Si=E_O_Si_ref.*T_ref./T;
    E_Si_Si=E_Si_Si_ref.*T_ref./T;
%     g_i_o=(4.29-16500/T);
    g_i_o=0;
   
KD_O=10.^(KD_O_p(1)+KD_O_p(2)./T+KD_O_p(3).*P./T+...
    +(g_i_o)./2.303...
    +E_O_O.*log(1-O_mol_frac_pick(1))./2.303...
    +E_O_Si./2.303.*Si_mol_frac_pick(1).*(1+log(1-Si_mol_frac_pick(1))./Si_mol_frac_pick(1)-1./(1-O_mol_frac_pick(1)))...
    -E_O_Si./2.303.*Si_mol_frac_pick(1).^2.*O_mol_frac_pick(1).*(1./(1-O_mol_frac_pick(1))+1./(1-Si_mol_frac_pick(1))+O_mol_frac_pick(1)./(2.*(1-O_mol_frac_pick(1)).^2)-1));

KD_Si=10.^(KD_Si_p(1)+KD_Si_p(2)./T+KD_Si_p(3).*P./T+E_Si_Si.*log(1-Si_mol_frac_pick(1))./2.303...
    +E_O_Si./2.303.*O_mol_frac_pick(1).*(1+log(1-O_mol_frac_pick(1))./O_mol_frac_pick(1)-1./(1-Si_mol_frac_pick(1)))...
    -E_O_Si./2.303.*O_mol_frac_pick(1).^2.*Si_mol_frac_pick(1).*(1./(1-Si_mol_frac_pick(1))+1./(1-O_mol_frac_pick(1))+Si_mol_frac_pick(1)./(2.*(1-Si_mol_frac_pick(1)).^2)-1));


KD_Ni=10.^(KD_Ni_p(1)+KD_Ni_p(2)./T+KD_Ni_p(3).*P./T);
% KD_O=10.^(KD_O_p(1)+KD_O_p(2)./T+KD_O_p(3).*P./T);
% 
% KD_Si=10.^(KD_Si_p(1)+KD_Si_p(2)./T+KD_Si_p(3).*P./T);

for w=1:numel(FeO_guess_v)
    
FeO_guess=FeO_guess_v(w);

NiO=FeO_guess.*(NiO_i+Ni_i)./((FeO_i+Fe_i-FeO_guess).*KD_Ni+FeO_guess);
Fe=FeO_i+Fe_i-FeO_guess;
Ni=NiO_i+Ni_i-NiO;
alpha=SiO2_i+Si_i;
gamma=Fe+Ni+FeO_i+NiO_i+3.*SiO2_i+O_i-FeO_guess-NiO+Si_i;
sigma=FeO_guess+NiO+MgO_i+AlO1_5_i+CaO_i;

quad_a=(3.*FeO_guess.^2-Fe.^2.*KD_Si);
quad_b=-(gamma.*FeO_guess.^2+3.*alpha*FeO_guess.^2+Fe.^2.*sigma.*KD_Si);
quad_c=alpha.*gamma.*FeO_guess.^2;
root_solve=roots([quad_a quad_b quad_c]);
SiO2=root_solve(2);
O=FeO_i+NiO_i+SiO2_i.*2+O_i-FeO_guess-NiO-2.*SiO2;
Si=SiO2_i+Si_i-SiO2;

%x=FeO;y=NiO,z=SiO2,u=MgO,m=Al1.5O,n=CaO
%a=Fe,b=Ni,c=O,d=Si

melt_mols=FeO_guess+NiO+SiO2+MgO_i+AlO1_5_i+CaO_i;
metal_mols=Fe+Ni+O+Si;

FeO_mol_frac(w)=FeO_guess./melt_mols;
NiO_mol_frac(w)=NiO./melt_mols;
SiO2_mol_frac(w)=SiO2./melt_mols;
MgO_mol_frac(w)=MgO_i./melt_mols;
AlO1_mol_frac(w)=AlO1_5_i./melt_mols;
CaO_mol_frac(w)=CaO_i/melt_mols;
Fe_mol_frac(w)=Fe./metal_mols;
Ni_mol_frac(w)=Ni./metal_mols;
O_mol_frac(w)=O./metal_mols;
Si_mol_frac(w)=Si./metal_mols;

melt_mols_logger=melt_mols;
metal_mols_logger=metal_mols;

KD_O_guess(w)=Fe_mol_frac(w).*O_mol_frac(w)./FeO_mol_frac(w);

iternum(w)=w;
KD_O_guess_logger(w)=KD_O_guess(w);


end

[min_comp,min_pos]=min(abs(KD_O_guess-KD_O));
min_pos_logger=min_pos;
min_comp_logger=min_comp;

FeO_mol_frac_pick=FeO_mol_frac(min_pos);
NiO_mol_frac_pick=NiO_mol_frac(min_pos);
SiO2_mol_frac_pick=SiO2_mol_frac(min_pos);
MgO_mol_frac_pick=MgO_mol_frac(min_pos);
AlO1_mol_frac_pick=AlO1_mol_frac(min_pos);
CaO_mol_frac_pick=CaO_mol_frac(min_pos);
Fe_mol_frac_pick=Fe_mol_frac(min_pos);
Ni_mol_frac_pick=Ni_mol_frac(min_pos);
O_mol_frac_pick=O_mol_frac(min_pos);
Si_mol_frac_pick=Si_mol_frac(min_pos);
KD_O_frac_er=(KD_O_guess(min_pos)-KD_O)./KD_O;
%FeO_guess_v=linspace(FeO_guess_v(min_pos)*0.95,FeO_guess_v(min_pos).*1.05,1e4);

FeO_guess_logger(:)=FeO_guess_v;
KD_O_logger=KD_O;

KD_Si_1_logger=KD_Si;
%KD_Si_2_logger=KD_Si_2;
KD_O_1_logger=KD_O;


%Values at the end
dIW=2.*log10(FeO_mol_frac_pick./Fe_mol_frac_pick);

FeO_guess_pass = FeO_guess_logger(:);

FeOmolfrac=FeO_mol_frac_pick;
NiOmolfrac=NiO_mol_frac_pick;
SiO2molfrac=SiO2_mol_frac_pick;
MgOmolfrac=MgO_mol_frac_pick;
AlO1molfrac=AlO1_mol_frac_pick;
CaOmolfrac=CaO_mol_frac_pick;
Femolfrac=Fe_mol_frac_pick;
Nimolfrac=Ni_mol_frac_pick;
Omolfracpick=O_mol_frac_pick;
Simolfrac=Si_mol_frac_pick;
% 
% 
A = [FeOmolfrac;NiOmolfrac;SiO2molfrac;Femolfrac;Nimolfrac;Omolfracpick; Simolfrac;...
    MgOmolfrac; AlO1molfrac; CaOmolfrac];
% 
% fontsize = 12;
% 
% 
% figure(1)
% subplot(3,3,1)
% plot(FeO_guess_v,SiO2_mol_frac,'-k',max(FeO_guess_v),melt_mol_frac_i(1),'ok', fontsize)
% title('SiO2', fontsize)
% subplot(3,3,2)
% plot(FeO_guess_v,FeO_mol_frac,'-k',max(FeO_guess_v),melt_mol_frac_i(3),'ok', fontsize)
% title('FeO', fontsize)
% subplot(3,3,3)
% plot(FeO_guess_v,NiO_mol_frac,'-k',max(FeO_guess_v),melt_mol_frac_i(6),'ok', fontsize)
% title('NiO')
% subplot(3,3,4)
% plot(FeO_guess_v,Fe_mol_frac,'-k',max(FeO_guess_v),metal_mol_frac_i(1),'ok', fontsize)
% title('Fe', fontsize)
% subplot(3,3,5)
% plot(FeO_guess_v,Ni_mol_frac,'-k',max(FeO_guess_v),metal_mol_frac_i(2),'ok', fontsize)
% title('Ni', fontsize)
% subplot(3,3,6)
% plot(FeO_guess_v,O_mol_frac,'-k',max(FeO_guess_v),metal_mol_frac_i(3),'ok', fontsize)
% title('O', fontsize)
% subplot(3,3,7)
% plot(FeO_guess_v,Si_mol_frac,'-k',max(FeO_guess_v),metal_mol_frac_i(4),'ok', fontsize)
% title('Si', fontsize)
