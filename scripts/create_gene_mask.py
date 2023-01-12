import os
import sys
import pandas as pd

refGene = int(os.environ['refGene'])
maskGenes = int(os.environ['maskGenes'])

print('Make sure all your files have ".txt" extension', maskGenes)
print('Bed-files will be generated here:', maskGenes)

