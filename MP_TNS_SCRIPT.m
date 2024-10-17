%************************************************************************************************************************
%***************************************COMPRESSION DE LA PAROLE PAR CODAGE EN SOUS-BANDES*******************************
%******************************************MP TNS Elaboré par Balti Chaima 1AMIndS2023***********************************
%************************************************************************************************************************

pkg load signal

%************************************************************************************************************************
%********************************************************Analyse*********************************************************
%************************************************************************************************************************

%************************************************Banc de filtres d’analyse***********************************************

%1)la lecture du signal audio en le stockant dans le signal s.
[s, Fe] = audioread("voix_homme_8.wav");
t=((0:size(s)(1)-1)*1/Fe)';

%2)filtre passe-bas d’ordre p = 7.
fc = (Fe/4);
p = 7;
HL=fir1(p,2*fc/Fe,'low',rectwin(p+1));
hb=impz(HL); %réponse impulsionnelle

%3)Représentation de la réponse fréquentielle (gain en amplitude et phase) du filtre passe-bas d’ordre p = 7.
figure(1)
[Hl,Wl]=freqz(HL,1);
subplot(2,1,1)
plot(Fe*Wl/(2*pi),abs(Hl),";Fenetre rectangulaire;");
title("Gain d'amplitude du filtre passe bas")
hold on
subplot(2,1,2)
plot(Fe*Wl/(2*pi),unwrap(angle((Hl)),pi),";Fenetre rectangulaire;");
hold on
title("Gain de phase du filtre passe bas")

%4)filtre passe-bas d’ordre p = 7 avec fenêtre de Hamming.
HLh=fir1(p,2*fc/Fe,'Low');
[Hlh,Wlh]=freqz(HLh,1);
subplot(2,1,1)
plot(Fe*Wlh/(2*pi),abs(Hlh),";Fenetre de Hamming;");
subplot(2,1,2)
plot(Fe*Wlh/(2*pi),unwrap(angle((Hlh)),pi),";Fenetre de Hamming;");

%5)filtre passe-bas d’ordre p = 24.
p=24;
HL24=fir1(p,2*fc/Fe,'Low');
figure(2)
[Hl24,Wl24]=freqz(HL24,1);
subplot(2,1,1)
plot(Fe*Wl24/(2*pi),abs(Hl24),";p=24;");
title("Gain d'amplitude du filtre passe bas pour différents ordres")
hold on
plot(Fe*Wl/(2*pi),abs(Hl),";p=7;");
subplot(2,1,2)
plot(Fe*Wl24/(2*pi),unwrap(angle((Hl24)),pi),";p=24;");
title("Gain de phase du filtre passe bas pour différents ordres")
hold on
plot(Fe*Wl/(2*pi),unwrap(angle((Hl)),pi),";p=7;");

%6)Application du filtre passe-bas sur le signal s et représentation du spectre du signal résultant, s1.
L=length(s);
s1=filter(HL24,1,s);
S1=fftshift(fft(s1,L));
freq=(-L/2 : (L-1)/2)*(Fe/L);
figure(3)
plot(freq,abs(S1),";Spectre du signal avec un filtre passe bas;")
hold on

%7)Représentation de la réponse fréquentielle (gain en amplitude et phase) du filtre passe-haut d’ordre p = 24.
k=0:p;
HH24=((-1).^k).*HL24;
figure(4)
[Hi,Wi]=freqz(HH24,1);
subplot(2,1,1)
plot(Fe*Wi/(2*pi),abs(Hi));
title("Gain d'amplitude du filtre passe haut")
hold on
subplot(2,1,2)
plot(Fe*Wi/(2*pi),unwrap(angle((Hi)),pi));
title("Gain de phase du filtre passe haut")
hold on

%8)Application du filtre passe-haut sur le signal s et représentation du spectre du signal résultant, s2.
s2=filter(HH24,1,s);
S2=fftshift(fft(s2,L));
figure(3)
plot(freq,abs(S2),";Spectre du signal avec un filtre passe haut;")

%9)signal reconstruit,s_r, à partir de s1 et s2.
s_r=s1+s2;
%soundsc(s);
%soundsc(s_r);

%10)Représentation dans le domaine temporel des signaux, s et sr.
figure(6)
subplot(2,1,1)

plot(s);
title("Le signal original dans le domaine temporelle")

subplot(2,1,2)
plot(s_r);
title("Le signal reconstitué dans le domaine temporelle")

%11)les filtres RII avec la fonction butter
[BH, AH] = butter (p, 2*fc/Fe,'high' );
[BL, AL] = butter (p, 2*fc/Fe,'low' );
s1i=filter(BL ,AL,s);
s2i=filter(BH ,AH,s);
s_ri=s1i+s2i;

figure(7)
subplot(2,1,1)
plot(s);
title("Le signal original dans le domaine temporelle")
subplot(2,1,2)
plot(s_ri);
title("Le signal reconstitué dans le domaine temporelle")


%******************************************************Décimation****************************************************

%Réduire la fréquence d’échantillonnage à celle de Nyquist
FeNyquist=Fe/2;

%1)Sous-échantillonnage des signaux s1 et s2
s1d = s1(1:2:end);
s2d = s2(1:2:end);

%2)Représentation dees spectres des signaux s1d et s2d.
Le=length(s1d);
freqe=(-Le/2 : (Le-1)/2)*(FeNyquist/Le);

S1d=fftshift(fft(s1d,Le));
figure(8)
plot(freqe,abs(S1d))
hold on
S2d=fftshift(fft(s2d,Le));
plot(freqe,abs(S2d),'r')

%*************************************************Quantification/Codage**********************************************

%1)fonction unifquant

%2)Application de unifquant sur les signaux s1d et s2d, pour différentes valeurs de l, choisies entre 2 et 12.
RSB1p=zeros(1,11);
RSB2p=zeros(1,12);
s1dq_vect=zeros(size(s1d)(1),11);
s2dq_vect=zeros(size(s2d)(1),11);

for l=2:12
  [s1dq_vect(:,l-1),RSB1p(1,l-1)]=unifquant(s1d,l);
  [s2dq_vect(:,l-1),RSB2p(1,l-1)]=unifquant(s2d,l);
endfor

KO=[2,4,6,12];
for l=1:4
  figure(9)
  subplot(2,2,l)
  plot(s1dq_vect(:,KO(l)-1))
  legend(['l = ' num2str(KO(l)-1) 'bits'])
  title("S1d quantifiée")
  hold on
  figure(10)
  subplot(2,2,l)
  plot(s2dq_vect(:,KO(l)-1))
  legend(['l = ' num2str(KO(l)-1) 'bits'])
  title("S2d quantifiée")
  hold
endfor

%choix de l1 et l2
%On prend l=7 pour les basses freq et 8 pour les hautes freq
l1=7;
l2=8;
L1=2^l1;
L2=2^l2;
[s1dq,RSB1]=unifquant(s1d,l1);
[s2dq,RSB2]=unifquant(s2d,l2);


%***************************************************Interpolation********************************************************

%1)suréchantillonnage des s1dq et s2dq
s1ech=zeros((size(s1dq)(1)*2)-1,1);
s2ech=zeros((size(s2dq)(1)*2)-1,1);

for i=1:2:(size(s1dq)(1)*2)-1
  s1ech(i,1)=s1dq((i+1)/2);
  s2ech(i,1)=s2dq((i+1)/2);
endfor

%2)Représentation des spectres des signaux résultants.
Le1=length(s1ech);
freq1=(-Le1/2 : (Le1-1)/2)*(Fe/Le1);
S1ech=fftshift(fft(s1ech,Le1));
figure(11)
plot(freq1,abs(S1ech))
hold on
S2ech=fftshift(fft(s2ech,Le1));
plot(freq1,abs(S2ech),'r')

%************************************************************************************************************************
%********************************************************Synthése********************************************************
%************************************************************************************************************************

%1)Application des filtres de synthèse sur s1ech et s2ech après suréchantillonnage et représentation de ses spectres.
s1echp=filter(HL24,1,s1ech);
s2echp=filter(HH24,1,s2ech);

S1echp=fftshift(fft(s1echp,Le1));
figure(12)
plot(freq1,abs(S1echp))
hold on
S2echp=fftshift(fft(s2echp,Le1));
plot(freq1,abs(S2echp),'r')

%2)Synthétisation du signal ss et calcul de la valeur du décalage temporel résultant.
ss=s1echp+s2echp;
% Calcul de la corrélation croisée entre s et ss
corr= xcorr(ss, s);
% Recherche de l'indice du maximum de la corrélation croisée
[max_corr, max_index] = max(corr);
%decalage temporel
decalage_temporel = (max_index-size(s)(1))/Fe;
%decalage temporel par bit
decalage_temporel_par_bit =decalage_temporel*(l1+l2)*Fe;
%Représentation
figure(13)
hold on
plot(t(1:1000)-decalage_temporel ,ss(1:1000),";signal ss;");
plot(t(1:1000),s(1:1000),";signal s;");
title('les signals ss et s dans le domaine temporelle');

%3)Pour différents choix de l1 et l2, procéder à l’écoute du signal synthétisé en le comparant au signal de départ, s.
D=zeros(11,11);
for l11=2:12
  for l22=2:12
    s1echl=zeros((size(s1dq_vect(l11-1))(1)*2)-1,1);
    s2echl=zeros((size(s2dq_vect(:,l22-1))(1)*2)-1,1);

    for i=1:2:(size(s1dq_vect(:,l11-1))(1)*2)-1
      s1echl(i,1)=s1dq_vect((i+1)/2,l11-1);
      s2echl(i,1)=s2dq_vect((i+1)/2,l22-1);
    endfor

    ssl=filter(HL24,1,s1echl)+filter(HH24,1,s2echl);
    %soundsc(s);
    %soundsc(ssl);
%4)Déterminer pour chaque choix, le débit binaire résultant de toute la chaîne ainsi que le RSB (deja calculer)
    %debit binaire D
    D(l11-1,l22-1)=Fe*(L1+L2);
  endfor
endfor

%5)l2=8 et l1=7
