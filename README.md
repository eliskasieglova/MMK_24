# MMK_24
Úloha do předmětu Matematická kartografie. Zpracovali Anna Brázdová, David Martínek a Eliška Sieglová.

Před spuštěním je nutné mít nainstalované knihovny _arcpy_, _numpy_, _pathlib_ , _geopandas_ a _os_. Také je nutné nastavit na řádku č. 7 _workspace_.

Pro generování glóbů využíváme dva scripty
- _map\_projections.py_: zde jsou ve slovníku uložené textové parametry gnomonické projekce pro každou plošku,
- _globes.py_: hlavní script, který generuje .pdf soubor s 12 ploškami dvanáctistěnového glóbu

Na začátku jsou 2 složky: _temp_ (zde jsou uložené dočasné soubory, složku je možné po proběhnutí skriptu smazat) a _results_ (zde jsou uložené výsledné shapefily, tuto složku nemazat).

Ve scriptu _globes.py_ jsou pro každou plošku vygenerované poledníky, rovnoběžky, a její okraj a ty jsou uložené jako shapefily do dočasné složky. Shapefily jsou převedeny do gnomonické projekce pomocí parametrů uložených ve scriptu _map\_projections.py_, protože každá ploška má jiný kartografický pól. _Reproject_ bylo provedeno pomocí knihovny _geopandas_, poněvadž knihovna _arcpy_ na _arcpy.Project\_management_ vracela _FatalError_. Výsledky reprojekce se ukládají do složky results jako shapefily pojmenované boundary/meridian/parallel_face{_číslo plošky_}.

Následně je vytvořen _Layout_, do kterého je vloženo 12 _map\_frame_, každý má tvar pravidelného pětiúhelníku. Do každého _map\_frame_ je vložena jedna stěna glóbu. Toto pdf je určeno k vytisknutí, vystřižení jednotlivých plošek, a slepení.




