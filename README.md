# MMK_24
Úloha do předmětu Matematická kartografie. Zpracovali Anna Brázdová, David Martínek a Eliška Sieglová.

Skript globes.py pomocí pro každou plošku dvanáctistěnu generuje její okraje a síť poledníků a rovnoběžek. Vytvoří 2 složky: temporary (zde jsou uložené dočasné soubory, složku je možné po proběhnutí skriptu smazat) a results (zde jsou uložené výsledné shapefily, tuto složku nemazat).

Před spuštěním je nutné mít nainstalované knihovny arcpy, numpy, pathlib a os. Také je nutné nastavit na řádku č. 7 workspace.

Výsledky se ukládají do složky results jako shapefily pojmenované boundary/meridian/parallel_face{_číslo plošky_}. 
