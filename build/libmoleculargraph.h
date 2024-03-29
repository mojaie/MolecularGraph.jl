char* smilestomol(char*, char*);
char* smartstomol(char*);
char* sdftomol(char*, char*);
int vertex_count(char*);
int edge_count(char*);
char* inchikey(char*);
double standard_weight(char*);
char* drawsvg(char*);
int drawpng(char*, char*, unsigned int, unsigned int);
int has_exact_match(char*, char*, char*);
int has_substruct_match(char*, char*, char*);
int tdmcis_size(char*, char*, char*);
int tdmces_size(char*, char*, char*);
double tdmcis_gls(char*, char*, char*);
double tdmces_gls(char*, char*, char*);

/* deprecated */
char* tdmces_gls_batch(char*);
char* smilestosvg(char*);
char* sdftosvg(char*);
double tdmcis_tanimoto(char*, char*, char*);
double tdmces_tanimoto(char*, char*, char*);
int tdmcis_dist(char*, char*, char*);
int tdmces_dist(char*, char*, char*);