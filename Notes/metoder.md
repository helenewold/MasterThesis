# Egenverdi-analyse
*Per 06.02.23*

1. Danne datasett med et sett punkter, 2 punkter med avstand d=10mm på dybder 20, 30 og 40 mm
2. Kjøre dette inn i getCapon
3. Hente ut R-matrisene fra getCapon dersom parametere skal regnes
4. Dermed regne eig-verdier og plotte alle disse langs en akse - midlertidig i midten mellom punktene.

Tanken videre er å dermed kjøre natt-kode, for å produsere resultater for færre punkter per sett, med 
forskjellige avstander mellom og dybder.

# Condition number
*Per 07.02.23*

Dette faller delvis under egenverdi-analyse. Innebærer stegene:
+ Plotte prosent av egenverdier under en viss grense - denne grensen er forskjellig fra simulert til ekte data. 
+ Kurve for å se på hvordan Condition Number endrer seg i et visst område. Dette faller i samme type uthenting av verdier som det egenverdi-analysen gjør.
Tanke videre: Vurdere om det trengs en viss ekstra robustifisering. Ikke gjort (7/2)