# CS 341 GTEx Project  

Team members: Nico Chaves, Noam Weinberger and Junjie (Jason) Zhu

## Structure

`/data`: includes metadata of processed data and example datasets; full datasets are stored on the server

`/preprocessing`: includes python scripts that used to process the expression data downloaded from [GTExPortal](http://www.gtexportal.org/home/datasets)

`/ipython_notebook`: includes ipython notebooks used to display main results of this project in an interactive fashion

## Instructions

`git pull`
`git add --all`
`git commit -m "MESSAGE"`
`git push`

## Data Preprocessing

We downloaded the Transcript RPKM file (`GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt.gz`) and meta-information (`GTEx_Data_V6_Annotations_SampleAttributesDS.txt`) from [GTExPortal](http://www.gtexportal.org/home/datasets) where the former includes expression values and the latter includes information about donor IDs and tissue types of each sample. Then we filtered out transcripts according to the following procedure:

1. Select transcripts that are mapped to genes in the GO database (list downloaded from [Ensembl](http://uswest.ensembl.org/biomart/martview/e9b91b8cc3de4a51e3a6f7cacad17699) Biomart)
2. Select top 10,000 transcripts with the highest variance across all samples in this dataset
3. 

|                     | Downloaded Transcript RPKM | After GO Term Filter | After Variance Filter |
| --------------------|:--------------------------:| :-------------------:|:---------------------:|
| Number of Variables | 195,747                    | 67,344               | 10,000                |


TODO: Write usage instructions

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## History

TODO: Write history

## Credits

TODO: Write credits

## License

TODO: Write license
