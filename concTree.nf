/*
 *   This file is part of 'concTree-NF'.
 *
 *   concTree-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   concTree-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with concTree-NF.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main concMSA-NF pipeline script
 *
 * @authors
 * Jia-Ming Chang <chang.jiaming@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Cedric Notredame <cedric.notredame@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 */


params.name          = "Concatenated MSAs and Consensus Phylogentic Trees Analysis"
params.input         = "$baseDir/tutorial/data/*.fa"
params.output        = "$baseDir/tutorial/results"
params.bootstraps    = 5
params.alignments    = 10
params.guidanceDir   = ''


log.info "c o n c T r e e  - N F  ~  version 0.1"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "input                  : ${params.input}"
log.info "output                 : ${params.output}"
log.info "alignments             : ${params.alignments}"
log.info "\n"


/*
 * Input parameters validation
 */

input_sequences         = file(params.input)
results_path            = file(params.output)

guidance_alignments_num = params.alignments as int


/*
 * Create a channel for input sequence files 
 */
 
datasets = Channel
    .fromPath( params.input )
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
    .map { file -> tuple( file.baseName, file ) }
 

process guidance2 {
    publishDir "$results_path/$datasetID/guidance2", mode: 'copy', overwrite: 'true'
  
    input:
    set val(datasetID), file(datasetFile) from datasets

    output:
    set val(datasetID), file ("selectedMSAs") into guidanceAlignments
    set val(datasetID), file ("base_alignment.aln") into baseAlignments   
    file("selectedMSAs") into alternativeAlignmentsDir
  
    script:
    //
    // Guidance2: Generating alternative alignments
    //
    
    """
    perl \$GUIDANCE2DIR/www/Guidance/guidance.pl \
        --seqFile ${datasetFile} \
        --msaProgram CLUSTALW \
        --seqType aa \
        --outDir \$PWD \
        --clustalw clustalw2

    cp MSA.CLUSTALW.aln.With_Names base_alignment.aln

    mkdir \$PWD/alternativeMSA
    tar -zxf MSA.CLUSTALW.Guidance2_AlternativeMSA.tar.gz -C alternativeMSA

    containsElement () {
        local e
        for e in "\${@:2}"; do [[ "\$e" == "\$1" ]] && return 0; done
        return 1
    }

    perl $baseDir/bin/hhsearch_cluster.pl alternativeMSA a3m $guidance_alignments_num

    while IFS=\\= read var value; do
        vars+=(\$var)
        values+=(\$value)
    done < repCluster.txt

    mkdir selectedMSAs

    ls alternativeMSA/ |while read file; 
    do
    if containsElement \$file \${vars[@]}
        then
            echo "Copying:  \$file"
            cp alternativeMSA/\$file selectedMSAs/.
        else
            echo "Skipping: \$file"
        fi
    done
    """
}

/*
 * Create a channel for the alternative alignments
 */

alternativeAlignmentsDir
    .flatMap { it.list() }
    .set { alternativeAlignmentFiles }

process trees {

    input:
    file(alternativeAlignment) from alternativeAlignmentFiles

    output:
    file("alternativeTree.nwk") into alternativeTreeFiles

    script:
    //
    // Generate Phylogenetic Trees: Generate phylogenetic trees in Newick format from each alternative alignment
    //

    """
    raxmlHPC-SSE3 -m PROTGAMMAJTT -p 4533 -T 1 -s ${basePhylip} -n tmp
    cp RAxML_bestTree.tmp alternativeTree.nwk
    rm *.tmp
    """
}
