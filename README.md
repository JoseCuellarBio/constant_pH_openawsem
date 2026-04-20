# constant_pH_openawsem

Implementacion experimental de un esquema de pH constante sobre OpenAWSEM. El flujo combina dinamica molecular en OpenMM/OpenAWSEM con intentos Monte Carlo de cambio de protonacion, actualizando en corrida el termino electrostatico de Debye-Huckel.

El directorio contiene tanto los scripts base como un caso de ejemplo ya preparado para la proteina `1ZUG`.

## Contenido del repositorio

### Raiz

- `mm_run.py`: script principal de simulacion.
- `forces_setup.py`: define el conjunto de fuerzas de OpenAWSEM y agrega el termino electrostatico dependiente de `charge.txt`.
- `pH_debyeHuckelTerms.py`: implementa el `CustomNonbondedForce` para Debye-Huckel con cargas por particula.
- `Montecarlo.py`: logica para seleccionar residuos ionizables, evaluar el entorno local y aceptar o rechazar cambios de protonacion.
- `mm_run_importaciones.py`: variante del runner principal.
- `helpers/generate_charge.py`: genera un archivo de cargas inicial a partir de una secuencia FASTA.
- `helpers/Graficos_pka_final.ipynb`: notebook de analisis.

### Ejemplo `1ZUG/`

El directorio [1ZUG](/home/matador/constant_pH_openawsem/1ZUG) incluye un sistema listo para correr y varios resultados ya generados:

- estructuras de entrada como `1ZUG-openmmawsem.pdb`, `crystal_structure.pdb` y `crystal_structure.fasta`
- parametros auxiliares de AWSEM como `single_frags.mem`, `single_frags.npy`, `gamma.dat`, `burial_gamma.dat`, `membrane_gamma.dat`, `ssweight`
- archivo de cargas inicial `charge.txt`
- copias locales de `mm_run.py`, `forces_setup.py`, `pH_debyeHuckelTerms.py` y `Montecarlo.py`
- script de analisis `mm_analyze.py`
- salidas de ejemplo como `movie.dcd`, `movie.pdb`, `output.log`, `info.dat`, `Hawsem.state`, `checkpnt.chk` y `time.dat`

## Como funciona

En `mm_run.py` hay dos modos principales:

- `-m 0`: dinamica a temperatura constante con interrupciones periodicas para intentar cambios de protonacion.
- `-m 1`: annealing de temperatura sin el ciclo de protonacion constante.

Cuando se usa `-m 0`, el flujo hace lo siguiente:

1. corre dinamica molecular durante `--interruptFrequency` pasos
2. construye un objeto `Proteina` con la geometria actual
3. elige un residuo ionizable al azar
4. evalua un cambio de carga con un criterio Monte Carlo que combina termino de pH y termino electrostatico local
5. si el cambio es aceptado, actualiza en contexto los parametros del ultimo force term, que corresponde a Debye-Huckel

El archivo `Hawsem.state` guarda el historial de estados de protonacion aceptados a lo largo de la corrida.

## Requisitos

Se asume un entorno de Python con:

- `openawsem`
- `openmm`
- `numpy`
- `pandas`
- dependencias habituales del stack de OpenAWSEM

Ejemplo de activacion:

```bash
conda activate openawsem
```

## Archivo de cargas

El termino electrostatico y el modo de pH constante esperan un archivo `charge.txt` en el directorio de trabajo. El formato es de dos columnas:

```text
resid charge
```

donde `resid` es el indice del residuo y `charge` es su carga inicial.

El helper [generate_charge.py](/home/matador/constant_pH_openawsem/helpers/generate_charge.py:1) permite construirlo desde una FASTA:

```bash
python helpers/generate_charge.py helpers/1ZUG.fasta -o charge.txt
```

Nota: el helper escribe `charges.txt` por defecto, mientras que la simulacion busca `charge.txt`.

## Ejecucion recomendada

Tal como esta organizado hoy el repositorio, la forma mas segura de correr el ejemplo es desde `1ZUG/`, porque alli estan juntos:

- los archivos de entrada de AWSEM
- `charge.txt`
- `forces_setup.py`
- `mm_analyze.py`

Ejemplo:

```bash
cd 1ZUG
python mm_run.py 1ZUG -s 10000 -r 100 --interruptFrequency 10 --pH 7.0 --tempStart 300 -m 0 -p CPU
```

Argumentos relevantes de `mm_run.py`:

- `protein`: nombre base del sistema
- `-s, --steps`: numero total de pasos
- `-r, --reportFrequency`: frecuencia de escritura de reportes
- `-m, --simulation_mode`: `0` para pH constante, `1` para annealing
- `--interruptFrequency`: pasos entre intentos de protonacion en modo `0`
- `--pH`: pH usado en el criterio Monte Carlo
- `--tempStart` y `--tempEnd`: temperaturas de inicio y fin
- `-p, --platform`: plataforma de OpenMM, por ejemplo `CPU`, `OpenCL` o `CUDA`
- `-t, --thread`: cantidad de threads si se usa `CPU`
- `-f, --forces`: archivo con la definicion de fuerzas
- `--fromCheckPoint`: reinicia desde un checkpoint

## Salidas tipicas

Una corrida puede producir:

- `output.log`: energia potencial y temperatura
- `movie.dcd`: trayectoria
- `movie.pdb`: snapshots en formato PDB
- `native.pdb`: estructura inicial reportada
- `checkpnt.chk`: checkpoint de OpenMM
- `time.dat`: tiempo total de simulacion
- `Hawsem.state`: historial de estados de protonacion
- `info.dat`: analisis energetico posterior sobre la trayectoria

## Limitaciones actuales

- El repositorio no esta empaquetado como modulo; depende de ejecutar los scripts desde un directorio con los archivos esperados por nombre.
- `mm_run.py` en la raiz invoca `mm_analyze.py` al finalizar, pero ese script no existe en la raiz del proyecto; la copia disponible esta dentro de `1ZUG/`.
- Hay duplicacion de scripts entre la raiz y `1ZUG/`. Si vas a mantener el proyecto en GitHub, conviene definir una sola fuente de verdad para evitar desfasajes.

## Archivos utiles para GitHub

Si queres publicar este proyecto tal como esta hoy, conviene revisar si queres versionar tambien los archivos pesados generados por una corrida, por ejemplo:

- `1ZUG/movie.dcd`
- `1ZUG/movie.pdb`
- `1ZUG/single_frags.npy`
- `1ZUG/checkpnt.chk`

Dependiendo de tu objetivo, puede tener mas sentido dejarlos fuera del repositorio y conservar solo inputs, scripts y una descripcion del flujo.
