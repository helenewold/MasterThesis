# Condition Number

Condition Number regnes nå ut i USTB postprocess getCapon, både før og
etter covarians matrisen regnes ut med diagonal loading. Condition Number
regnes vha. matlab-funksjonen cond().

Det tas (ikke, det ble tidligere tatt) høyde for at condition nummeret blir uendelig stort. I de
tilfellene settes verdien til 0, for å unngå tilfeller av NaN eller ikke
plottbare resultater.

Når condition number er utregnet, kan det hentes ut slik:
'''
% Initier uff-objektet
b_data_CondNr                   =   uff.beamformed_data(b_data_mv_getCapon);
b_data_CondNr_after             =   uff.beamformed_data(b_data_mv_getCapon);
% Hent ut dataen
b_data_CondNr.data(:,:)         =   mv_getCapon.CN;
b_data_CondNr_after.data(:,:)   =   mv_getCapon.CN_after;
'''

Dette er da med forbehold om at det eksisterer et beamforming-objekt med
navn b_data_mv_getCapon, og at getCapon structuren kalles mv_getCapon.


Ved bruk av USTB for å plotte, ved
'''
fig2 = figure(2);
b_data_CondNr.plot(fig2, 'Condition number - Before diagonal loading', 400);

fig3 = figure(3);
b_data_CondNr_after.plot(fig3, 'Condition number - After diagonal loading', 400);
'''

får vi følgende resultat. Dette resultatet er altså plottet med logaritmisk
akse.

<p float="left">
  <img src=../Figures/Condition_Number/USTB_a_first.png width=40% />
  <img src=../Figures/Condition_Number/USTB_b_first.png width=40% />
</p>

Videre kan getScanConvertedImage brukes, ikke med logaritmisk akse.

<p float="left">
  <img src=../Figures/Condition_Number/getScan_first.png width=75% />
</p>

Resultatene er i skrivende stund (10.11.22) udiskuterbare, da enkelte
verdier er uendelig store, eller har størrelsesorden opp til 1e15.

Hva anses som et for stort condition number?

## Edit per 16.11:
Condition number plottes vha. getScan. Først deles condition numrene inn i
seksjoner mht størrelsen på verdien. Resultatet vises her:

<p float="left">
  <img src=../Figures/Condition_Number/CondNr_Sections.png width=75% />
</p>

Her er fargene definert i 0-2 altså bakgrunnen til scan-plottet.
2-4  representerer verdier mindre enn 100. 
4-6  representerer verdier mellom 100 og 10 000 000.
6-8  representerer verdier mellom 10 000 000 og 10 000 000 000 000 000.
8-10 representerer verdier større enn 10 000 000 000 000, altså 
hovedsakelig uendelige verdier.


Her ser vi altså at en seksjon mellom 5-15 mm dybde er uendelige verdier,
hovedsakelig de tilfellene der kovariansmatrisen ikke er invertibel. Dette
er et resultat av at dataen er simulert, som i enkelte tilfeller kan gi 0
på en pixelverdi. Da vil altså kovariansmatrisen bli en null-matrise, som
ikke er invertibel.






