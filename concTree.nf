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
params.bootstraps    = 100
params.alignments    = 100
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
bootstraps = params.bootstraps as int

/*
 * Create a channel for input sequence files 
 */
 
datasets = Channel
    .fromPath( params.input )
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
    .map { file -> tuple( file.baseName, file ) }
 

process guidance2 {
    tag "guidance2: $datasetID"

    publishDir "$results_path/$datasetID/guidance2", mode: 'copy', overwrite: 'true'
  
    input:
    set val(datasetID), file(datasetFile) from datasets

    output:
    set val(datasetID), file ("selectedMSAs") into guidanceAlignments
    set val(datasetID), file ("base_alignment.aln") into baseAlignments, baseAlignmentsForBootstrap   
    set val(datasetID), file("selectedMSAs/*") into alternativeAlignments mode flatten
  
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

process concatenate {
    publishDir "$results_path/$datasetID/alignments", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), file(alternativeAlignmentsDir) from guidanceAlignments
    set val(datasetID), file(base_alignment) from baseAlignments

    output:
    set val(datasetID), file("base_alignment.phylip") into baseAlignmentPhylips
    set val(datasetID), file("concatenated_alignment.phylip") into concatenatedAlignmentPhylips
    set val(datasetID), file("base_alignment_concatenated.phylip") into baseAlignmentConcatenated

    script:
    //
    // Concatenate Process: Generate the base and concatenated alignment in PHYLIP format
    // Also create a control alignment containing the base alignment concatenated 
    // 

    """
    shopt -s nullglob
    alternativeMSAsFASTA=( ${alternativeAlignmentsDir}/*".fasta" )
    alternativeMSAsPHYLIP=("\${alternativeMSAsFASTA[@]/.fasta/.phylip}")

    x=0
    mkdir baseAlignmentConcatenated
    for i in "\${alternativeMSAsFASTA[@]}"
    do
        outfileBase=\${i%.fasta}
        esl-reformat phylip \$i > \$outfileBase.phylip
        esl-reformat afa ${base_alignment} > baseAlignmentConcatenated/baseAlignment_\$x.fa
        let x=x+1
    done

    baseMSAsFASTA=( "baseAlignmentConcatenated/baseAlignment_*.fa" )
    concatenate.pl --aln \${baseMSAsFASTA[@]} --out base_alignment_concatenated.aln
    esl-reformat phylip base_alignment_concatenated.aln > base_alignment_concatenated.phylip
  
    concatenate.pl --aln \${alternativeMSAsFASTA[@]} --out concatenated_alignment.aln
    esl-reformat phylip concatenated_alignment.aln > concatenated_alignment.phylip
    esl-reformat phylip ${base_alignment} > base_alignment.phylip
    """
}

process trees {
    tag "trees: $datasetID"
    publishDir "$results_path/$datasetID/trees", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), file(basePhylip) from baseAlignmentPhylips
    set val(datasetID), file(concatenatedPhylip) from concatenatedAlignmentPhylips
    set val(datasetID), file(baseAlignmentConcatenatedPhylip) from baseAlignmentConcatenated

    output:
    set val(datasetID), file("base_tree.nwk") into baseTrees
    set val(datasetID), file("conc_tree.nwk") into concatenatedTrees
    set val(datasetID), file("base_conc_tree.nwk") into baseConcatenatedTrees

    script:
    //
    // Generate Phylogenetic Trees: Generate the base and concatenated phylogenetic trees in Newick format
    //

    """
    raxmlHPC-SSE3 -m PROTGAMMAJTT -p 4533 -T 1 -s ${basePhylip} -n tmp
    cp RAxML_bestTree.tmp base_tree.nwk
    rm *.tmp
    
    raxmlHPC-SSE3 -m PROTGAMMAJTT -p 4533 -T 1 -s ${concatenatedPhylip} -n tmp
    cp RAxML_bestTree.tmp conc_tree.nwk
    rm *.tmp

    raxmlHPC-SSE3 -m PROTGAMMAJTT -p 4533 -T 1 -s ${baseAlignmentConcatenatedPhylip} -n tmp
    cp RAxML_bestTree.tmp base_conc_tree.nwk
    rm *.tmp
    """
}

process paramastrap_trees {
    tag "paramastrap trees: $datasetID"

    input:
    set val (datasetID), file(alternativeAlignment) from alternativeAlignments

    output:
    set val (datasetID), file("paramastrap.phylip") into paramastrapPhylips
    set val (datasetID), file("paramastrap.nwk") into paramastrapTrees

    script:
    //
    // Generate Paramastrap Trees: Generate paramastrap trees in Newick format from each alternative alignment
    //

    """
    seed=4533
    esl-reformat phylip ${alternativeAlignment} > alternativeAlignment.phylip
    echo -e "alternativeAlignment.phylip\nR\n1\nY\n\$seed\n" | seqboot
    mv outfile paramastrap.phylip

    raxmlHPC-SSE3 -m PROTGAMMAJTT -p 4533 -T 1 -s paramastrap.phylip -n tmp
    cp RAxML_bestTree.tmp paramastrap.nwk
    rm *.tmp
    """
}

process bootstrap_trees {
    tag "bootstrap trees: $datasetID"

    input:
    each x from 1..bootstraps
    set val (datasetID), file(baseAlignment) from baseAlignmentsForBootstrap

    output:
    set val (datasetID), file("bootstrap.phylip") into bootstrapPhylips
    set val (datasetID), file("bootstrap.nwk") into bootstrapTrees

    script:
    //
    // Generate Bootstrap Trees: Generate bootstrap trees in Newick format from the base alignment
    //

    """
    seed=4533
    esl-reformat phylip ${baseAlignment} > baseAlignment.phylip
    echo -e "baseAlignment.phylip\nR\n1\nY\n\$seed\n" | seqboot
    mv outfile bootstrap.phylip

    raxmlHPC-SSE3 -m PROTGAMMAJTT -p 4533 -T 1 -s bootstrap.phylip -n tmp
    cp RAxML_bestTree.tmp bootstrap.nwk
    rm *.tmp
    """
}


paramastrapTrees
    .collectFile() { item ->
       [ "${item[0]}.paramastraps.txt", item[1] ]
    }
    .set { catParamastrapTrees }


bootstrapTrees
    .collectFile() { item ->
       [ "${item[0]}.bootstraps.txt", item[1] ]
     }
    .set { catBootstrapTrees }

process nodeSupport {
    publishDir "$results_path/$datasetID/nodeSupport/", mode: 'copy', overwrite: 'true'

    input:
    file(catBootstrap) from catBootstrapTrees
    file(catParamastrap) from catParamastrapTrees
    set val(datasetID), file(base_tree) from baseTrees
    set val(datasetID), file(conc_tree) from concatenatedTrees
    set val(datasetID), file(base_conc_tree) from baseConcatenatedTrees

    output:
    set val(datasetID), file('nodeSupportForBaseTree.result') into nodeSupportForBaseTrees
    set val(datasetID), file('nodeSupportForConcTree.result') into nodeSupportForConcTrees
    set val(datasetID), file('nodeSupportForBaseConcTree.result') into nodeSupportForBaseConcTrees

    script:
    """ 
    if [[ $datasetID =~ ^tips([0-9]+)_([0-9].[0-9])_[0-9]+.[0-9]+\$ ]]
    then
        reference_tree=${baseDir}/tutorial/data/ref_trees/tips\${BASH_REMATCH[1]}/asymmetric_\${BASH_REMATCH[2]}.unroot.tree
    else
        echo "unable to parse string ${datasetID}"
    fi 

    echo ">${catBootstrap}" > bootstrapReplicateFileList.txt
    echo ">${catParamastrap}" >> bootstrapReplicateFileList.txt

    t_coffee -other_pg seq_reformat -in ${base_tree} -in2 bootstrapReplicateFileList.txt -action +tree2ns \$reference_tree > nodeSupportForBaseTree.result

    t_coffee -other_pg seq_reformat -in ${conc_tree} -in2 bootstrapReplicateFileList.txt -action +tree2ns \$reference_tree > nodeSupportForConcTree.result

    t_coffee -other_pg seq_reformat -in ${base_conc_tree} -in2 bootstrapReplicateFileList.txt -action +tree2ns \$reference_tree > nodeSupportForBaseConcTree.result

    """
}
