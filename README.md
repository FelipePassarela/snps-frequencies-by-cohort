# Frequências genotípicas por coorte

Pacote Python gerenciado por uv, com código em `src/usp_gene_correlation`.
Calcula frequências de um VCF e acrescenta `name` e `rsid` por posição.

## Instalação e execução

Python 3.13+. Use Linux, macOS ou WSL para a dependência pysam.
No Windows, execute os comandos dentro do WSL; não compartilhe um ambiente
virtual Windows com o ambiente Linux.

```sh
uv sync
uv run usp-gene-correlation data/longa.vcf "data/cohorts/WES CEGH x UFES(INFECTADOS).csv" -o output/longa
uv run usp-gene-correlation data/grave.vcf "data/cohorts/WES CEGH x UFES(INFECTADOS).csv" -o output/grave
uv run pytest
```

Aceita vários arquivos de coorte. `--output-dir`, `--output_dir` e `-o` são equivalentes;
sem saída, imprime o resultado. Também é possível usar
`uv run python -m usp_gene_correlation ...`.

## Anotações

O catálogo padrão está em `src/usp_gene_correlation/data/snps_annotations.csv` e é
incluído no pacote instalado. Contém 44 SNPs: 19 transcritos da imagem
fornecida pelo usuário (longa) e 25 da tabela enviada na conversa (grave).
Os nomes foram preservados como fornecidos, sem validação em API externa.

Para substituir o catálogo:

```sh
uv run usp-gene-correlation data/longa.vcf coorte.csv --annotations meus_snps.csv -o output
```

Formato obrigatório:

```csv
chr_pos,name,rsid
chr10_69338539,HK1,rs7912524
```

As posições são de base 1 e precisam corresponder à mesma montagem do VCF.
Os arquivos locais têm referências a hg38 no cabeçalho; o catálogo fornecido
deve ser usado com esses dados. Não há conversão entre montagens nem
desambiguação por alelos. `10_69338539` e `chr10_69338539` são equivalentes.
Posições duplicadas, campos incompletos e rsIDs inválidos são rejeitados.
Variantes não encontradas são preservadas com campos vazios e aviso.

`SnpAnnotationSource` define `lookup(chr_pos) -> SnpAnnotation | None`.
A implementação `CsvAnnotationSource` carrega e valida o catálogo CSV.
Um adaptador de API pode implementar esse mesmo contrato e ser injetado em
`main(argv, annotation_source=fonte)` ou `annotate_frequencies(df, fonte)`;
nenhuma consulta de rede acontece durante a execução padrão.

## Entradas e saída

Os CSVs de coorte podem ter cabeçalho; a coluna de IDs é detectada pelo padrão
`C01234-ExC123-xgenV1`. Cabeçalhos e outros valores são descartados, assim
como IDs duplicados. Somente amostras presentes no VCF são incluídas.

A saída preserva a ordem das variantes e as colunas de frequência existentes:

```text
variant,name,rsid,0/0,0/1,1/1,./.
```

O denominador inclui todas as amostras selecionadas, inclusive chamadas
ausentes. Genótipos parcialmente ausentes são contados como `./.`;
heterozigotos invertidos e faseados são normalizados. Categorias ausentes
recebem zero. Genótipos chamados não diploides ou multialélicos são
rejeitados para evitar perder frequências silenciosamente. Coortes sem
amostras correspondentes ou VCFs sem variantes geram erro.

## Organização

- `cli.py`: CLI Typer com validação de caminhos e execução por coorte.
- `cohorts.py`: objeto imutável `Cohort(name, sample_ids)` e leitura de CSVs
  por `read_cohorts(paths) -> list[Cohort]`. Os IDs são tuplas, sem preenchimento
  de linhas ausentes entre coortes de tamanhos diferentes.
- `genotypes.py`: leitura VCF e frequências.
- `annotations/`: pacote com API pública em `__init__.py`, modelo e normalização em
  `models.py`, contrato em `base.py`, fonte CSV e catálogo padrão em `csv.py` e enriquecimento em `enrichment.py`.
- `src/usp_gene_correlation/data/snps_annotations.csv`: catálogo distribuído.
- `tests/conftest.py`: fixtures compartilhadas de VCF, coorte e catálogo sintéticos.
- `tests/test_annotations.py`, `test_cohorts.py`, `test_genotypes.py` e
  `test_cli.py`: testes separados por responsabilidade, incluindo integração da CLI.
- `scripts/verify_local_outputs.py`: verificação manual dos dados locais;
  atualiza os CSVs apenas após confirmar que as frequências foram preservadas.

O `uv.lock` fixa as dependências; `uv sync --locked` reproduz o ambiente.
`requirements.txt` é uma exportação para ferramentas que usam pip.
Para atualizá-lo: `uv export --no-dev --no-emit-project --format requirements-txt --output-file requirements.txt`.
