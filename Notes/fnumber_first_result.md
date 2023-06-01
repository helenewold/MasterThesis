# Results from fnumber and apodization matrix implementation in getCapon
Endringene i getCapon er å legge til 
```
% If apodization matrix is given
rx_apod = squeeze(apod_matrix(i, j, :));
idx = find( abs(rx_apod) > 0.16 );
M_new = length(idx);
L_new = floor(L_frac*M_new);
I_new = eye(L_new);
a = ones(L_new, 1);
```
i den doble for-loopen for å anvende apodiseringsmatrisen. Gjennom koden
etter dette er altså f.eks. L byttet ut med L_new, for å oppdatere med nye tall.
Dette gjennomføres kun når apod_matrix eksisterer.

Resultatene blir som følger: 
<p float="left">
  <img src=../Figures/Capon_Section/fnumber_getCapon.png width=40% />
  <img src=../Figures/Capon_Section/fnumber_USTB.png width=40% /> 
</p>
Ved analyse ser man hvordan resultatet er visuelt likt, annet enn i kantene
der USTB tar hensyn til nTimeAverage, mens getCapon ignorerer.

Ganger med length(idx) på slutten - hvorfor?
