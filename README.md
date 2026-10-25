# Frequências genotípicas por coorte

Calcula frequências de genótipos a partir de um VCF, separadas por coorte,
e acrescenta o nome e o rsID de cada SNP.

## Instalação e uso

Requer Python 3.13+ e uv. Use Linux, macOS ou WSL por causa da dependência pysam.

```sh
uv sync --locked
uv run usp-gene-correlation data/longa.vcf "data/cohorts/WES CEGH x UFES(INFECTADOS).csv" -o output/longa
```

Aceita vários CSVs de coorte após o caminho do VCF. Gera um CSV por coorte,
usando o nome do arquivo de entrada. `--output-dir`, `--output_dir` e `-o`
são equivalentes. Sem essa opção, imprime a tabela em stdout; mensagens
de progresso e avisos vão para stderr.

Para consultar todas as opções:

```sh
uv run usp-gene-correlation --help
```

## Entradas

O VCF deve conter os genótipos e os nomes das amostras. Cada CSV de coorte
deve ter uma coluna com IDs como `C01234-ExC123-xgenV1`. Cabeçalhos são
opcionais; valores fora desse padrão e IDs duplicados são descartados.
Somente amostras presentes no VCF entram no cálculo.

O catálogo incluído contém 44 SNPs dos conjuntos longa e grave.
Para usar outro catálogo:

```sh
uv run usp-gene-correlation dados.vcf coorte.csv --annotations meus_snps.csv -o output
```

O CSV de anotações deve seguir este formato:

```csv
chr_pos,name,rsid
chr10_69338539,HK1,rs7912524
```

As posições são de base 1 e devem usar a mesma montagem do VCF. O catálogo
incluído foi fornecido para os dados locais em hg38, sem validação externa.
Não há conversão entre montagens nem distinção por alelos na mesma posição.
`10_69338539` e `chr10_69338539` são equivalentes. Posições duplicadas,
campos incompletos e rsIDs inválidos são rejeitados. Variantes sem anotação
permanecem na saída com nome e rsID vazios, acompanhadas de um aviso.

## Resultados

As linhas seguem a ordem das variantes no VCF:

```text
variant,name,rsid,0/0,0/1,1/1,./.
```

As frequências incluem todas as amostras selecionadas no denominador,
inclusive genótipos ausentes, e somam 1 por variante.
`0/0` representa homozigotos de referência, `0/1` heterozigotos,
`1/1` homozigotos alternativos e `./.` chamadas ausentes.

Chamadas parcialmente ausentes são contadas como `./.`; heterozigotos
invertidos e faseados são normalizados. Categorias sem ocorrências recebem
zero. Chamadas não diploides ou multialélicas são rejeitadas. Coortes sem
amostras correspondentes ou VCFs sem variantes geram erro.

## Testes

```sh
uv run pytest
```

Diagrama da arquitetura: [PlantUML](docs/architecture.puml).
