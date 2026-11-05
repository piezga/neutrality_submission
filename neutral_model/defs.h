


typedef struct {
  int ix,iy;         // posizione attuale
  int name;           // indicizza il walker con la condizione iniziale
  int nameptr;        // indicizza con la condizione iniziale del walker con cui ha colliso
  int specie;         // indice della specie
} walker_type;

typedef struct {
  int nameptr;    // indicizza con la condizione iniziale del walker con cui ha colliso
  int specie;     // indice della specie -1 se uno di quelli eliminati
} species_type;

typedef struct {
  int size,numdata;
  int *list;    
} table_type;
