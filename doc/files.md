# DiGeST : Files

##global.R

Global.R is the first file executed at the launch of the Shiny applcation. It contains initialization and generic functions. In particular:

* It is where the environment is set for a cluster or a standalone use
* The CliniPhenome section has methods to update the phenotype DB from CliniPhenome, and to query the phenotype DB
* The Highlander section allows to determine what columns are ket from Highlander
* The generic functions are for displaying a table (getWidgetTable), loading variant data (loadData) and loading a set of results (procRes)

##UI.R

The user interface. 

* It is made of four tabs: Phenotype manager, Gene & variant filtering manager, Scoring tool, and Results explorer. 
* Filtering components are in the filterPhenotypes and filterVarints file (see below). 
* Result components are in the server file (see below). 
* The UI also has a small javascript component to deal with the GET request (see API doc for details).

##server.R

The server functions. Initialization of global variables:

* sessionvalues is used to store all Shiny reactive values
* data contains variant data. Initialized with empty filter
* nbRowsExceededWarningMessage: used to emit warning if more than 1000 variants retrieved
* phenotypes contain phenotype data. Initialized with empty filter
* analysesNames contains the set of analyses for a user
* variantDataGene: Used in results. Contains the set of variants in a given gene
* logged_user: The logged user


The server file is divided in the following parts 

* User login and UI controls
* Phenotypes group manager
* Variant group manager
* Ranking engine
* Results explorer

##filterPhenotypes.R and filterVariants.R

Contain the filtering widget logic for phenotypes and variants. Both file have the (almost) same structure, allowing to create, load, save, delete and apply filters. Function names are adapted so phenotype and variant filters can be differentiated. The variant filter further contain a function to be initialized with the sample_id GET call (see API).
