# enrichment
Simple script to perform set enrichment analysis

`ses.py` allows to perform a gene set enrichment analysis from the command line. It takes as input the gene set and will compare it to all the sets in a file containing predefined sets (usually a set corresponds to genes sharing an annotation, e.g. the genes annotated with the [Gene Ontology](http://geneontology.org) term GO:0006418 *tRNA aminoacylation for protein translation*	from the Biological Process branch).

## set comparisons

The **query set** is compared to each **target set** through a **binomial test**, thus each comparison yields a *p-value* which can be interpreted as the probability to obtain by chance at least the number of elements in common between the two compared sets (considered as two random samples from the same population).

# Example
The following command line will consider the set of identifiers {ALAS, ARGS, ASNS, ASPS, CYSS, GLTX, GLYQ, GLYS, HISS, ILES} and search for most *similar* sets in the file `data/EcolA.go.set`. The *p-values* are adjust by FDR because of the multiple testing performed and the theshold for *p-values* is set to 0.05.

```sh
./ses.py --query 'ALAS ARGS ASNS ASPS CYSS GLTX GLYQ GLYS HISS ILES' --sets data/EcolA.go.sets --adjust --alpha 0.05
```

The begining of the output produced is:
<pre>
GO:0006418	3.3913752513694115e-22	10/26	biological_process: tRNA aminoacylation for protein translation	GLTX, GLYQ, GLYS, ARGS, ASNS, ILES, ALAS, ASPS, HISS, CYSS
GO:0043039	4.946294932473471e-22	10/27	biological_process: tRNA aminoacylation	GLTX, GLYQ, GLYS, ARGS, ASNS, ILES, ALAS, ASPS, HISS, CYSS
GO:0016875	4.946294932473471e-22	10/27	molecular_function: ligase activity, forming carbon-oxygen bonds	GLTX, GLYQ, GLYS, ARGS, ASNS, ILES, ALAS, ASPS, HISS, CYSS
GO:0016876	4.946294932473471e-22	10/27	molecular_function: ligase activity, forming aminoacyl-tRNA and related compounds	GLTX, GLYQ, GLYS, ARGS, ASNS, ILES, ALAS, ASPS, HISS, CYSS
GO:0004812	4.946294932473471e-22	10/27	molecular_function: aminoacyl-tRNA ligase activity	GLTX, GLYQ, GLYS, ARGS, ASNS, ILES, ALAS, ASPS, HISS, CYSS
GO:0043038	7.115782749919615e-22	10/28	biological_process: amino acid activation	GLTX, GLYQ, GLYS, ARGS, ASNS, ILES, ALAS, ASPS, HISS, CYSS
GO:0006399	5.3165129894909594e-17	10/86	biological_process: tRNA metabolic process	GLTX, GLYQ, GLYS, ARGS, ASNS, ILES, ALAS, ASPS, HISS, CYSS
GO:0016874	8.376593900782695e-17	10/90	molecular_function: ligase activity	GLTX, GLYQ, GLYS, ARGS, ASNS, ILES, ALAS, ASPS, HISS, CYSS
GO:0034660	1.2573681058420214e-15	10/118	biological_process: ncRNA metabolic process	GLTX, GLYQ, GLYS, ARGS, ASNS, ILES, ALAS, ASPS, HISS, CYSS
...
</pre>
