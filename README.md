# MMK_24
Úloha do předmětu Matematická kartografie. Zpracovali Anna Brázdová, David Martínek a Eliška Sieglová.

Skript _globes.py_ pomocí pro každou plošku dvanáctistěnu generuje její okraje a síť poledníků a rovnoběžek. Vytvoří 2 složky: _temporary_ (zde jsou uložené dočasné soubory, složku je možné po proběhnutí skriptu smazat) a _results_ (zde jsou uložené výsledné shapefily, tuto složku nemazat).

Před spuštěním je nutné mít nainstalované knihovny _arcpy_, _numpy_, _pathlib_ a _os_. Také je nutné nastavit na řádku č. 7 _workspace_.

Výsledky se ukládají do složky results jako shapefily pojmenované boundary/meridian/parallel_face{_číslo plošky_}. Ty je možné zobrazit v ArcGISu. Pro zobrazení výsledků je nutné mít pro každou plošku otevřenou jinou _Map_, které se nastaví odpovídající projekce North/South gnomonic a s odpovídajícím kartografickým pólem dle následující tabulky.


| Face | u | v | 
| ------------- | ------------- | ------------- |
| 1	| -90 |	0 |
| 2	| 90	| 0 |
| 3	| -20.90515	| 36 |
| 4	| -20.90515	| 108 |
| 5	| -20.90515	| 180 | 
| 6	| -20.90515	| 252 |
| 7	| -20.90515	| 324 |
| 8	| 31.71745	| 0 |
| 9	| 31.71745	| 72 |
| 10	| 31.71745	| 144 |
| 11	| 31.71745	| 216 |
| 12	| 31.71745	| 288 |


Náš výsledek připravený k tisku je uložen v této složce pod jménem _globe_faces.png_.
