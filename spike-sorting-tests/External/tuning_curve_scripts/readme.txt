this is code from marcel stimberg, mail from 22.06.2012

Hi Felix,

anbei scripts, die Gauß oder von Mises Funktionen fitten, mit (_fixed_dist) oder ohne (_free) Annahme, dass die Peaks 180° entfernt sind. Die fit-Funktionen erwarten eine Matrix mit den Tuning-Kurven (eine Zelle pro Zeile) -- im Skript steht die Auflösung der Tuningkurven (22.5°) fest drin, das musst Du evtl. noch anpassen. Ach so, in der _chisquare-Variante wird noch die Genauigkeit der einzelnen Datenpunkte berücksichtigt, was man eigentlich immer machen sollte, ich aber manchmal mangels dafür nötiger Daten nicht machen konnte...

Viel Spaß damit,
  Marcel