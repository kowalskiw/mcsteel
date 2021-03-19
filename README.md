# McSteel
Set of scripts for stochastic analysis of structure exposed to a fire. As a submodel McSteel uses SAFIR® software by Univerity of Liège (https://www.uee.uliege.be/cms/c_4016386/en/safir). Multisimulation approach is used (https://bibliotekasgsp.locloudhosting.net/items/show/120).


* `mc.py` - generating multisimulation scenarios (Safir Thermal 2D);

* `multi.py` - calculating already generated scenarios (Safir Thermal 2D);

* `fires.py` - fire models database;

* `chart.py` - drawing results distributions;

* `export.py` - saving results to the database;

* `fdsafir.py` - transforming ISO fire load to LOCAFI and running single simulation (Safir Structural 3D & Safir Thermal 2D);

* `selsim.py` - selecting the severest percentile of scenarios, preparing and running global analysis of structure for those (Safir Structural 3D).
