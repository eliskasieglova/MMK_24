# MMK_24
Úloha do předmětu Matematická kartografie. Zpracovali Anna Brázdová, David Martínek a Eliška Sieglová.

Skript _globes.py_ pomocí pro každou plošku dvanáctistěnu generuje její okraje a síť poledníků a rovnoběžek. Vytvoří 2 složky: _temporary_ (zde jsou uložené dočasné soubory, složku je možné po proběhnutí skriptu smazat) a _results_ (zde jsou uložené výsledné shapefily, tuto složku nemazat).

Před spuštěním je nutné mít nainstalované knihovny _arcpy_, _numpy_, _pathlib_ a _os_. Také je nutné nastavit na řádku č. 7 _workspace_.

Výsledky se ukládají do složky results jako shapefily pojmenované boundary/meridian/parallel_face{_číslo plošky_}. Ty je možné zobrazit v ArcGISu. Pro zobrazení výsledků je nutné mít pro každou plošku otevřenou jinou _Map_, které se nastaví odpovídající projekce North/South gnomonic a s odpovídajícím kartografickým pólem dle následující tabulky.






Náš výsledek připravený k tisku je uložen v této složce pod jménem _globe_faces.png_.
