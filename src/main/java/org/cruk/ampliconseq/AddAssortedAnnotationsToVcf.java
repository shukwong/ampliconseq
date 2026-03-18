/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.ampliconseq;

import java.io.File;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Utility for adding various annotations to variants in a VCF file.
 *
 * @author eldrid01
 */
@Command(name = "add-assorted-annotations-to-vcf", versionProvider = AddAssortedAnnotationsToVcf.class, description = "\nAdds various annotations to variants in a VCF file.\n", mixinStandardHelpOptions = true)
public class AddAssortedAnnotationsToVcf extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger();

    // annotation names
    private static final String MULTIALLELIC = "MULTIALLELIC";
    private static final String FIVE_PRIME_SEQUENCE_CONTEXT = "FivePrimeContext";
    private static final String THREE_PRIME_SEQUENCE_CONTEXT = "ThreePrimeContext";
    private static final String INDEL_LENGTH = "IndelLength";
    private static final String HOMOPOLYMER_LENGTH = "HomopolymerLength";

    @Option(names = { "-i", "--input" }, required = true, description = "Input VCF file (required).")
    private File inputVcfFile;

    @Option(names = { "-r",
            "--reference-sequence" }, required = true, description = "Reference sequence FASTA file which must be indexed and have an accompanying dictionary (required).")
    private File referenceSequenceFile;

    @Option(names = { "-o",
            "--output" }, required = true, description = "Output VCF file with annotated varaints (required).")
    private File outputVcfFile;

    @Option(names = {
            "--sequence-context-length" }, description = "The number of bases of sequence context to record for both 5' and 3' context (default: ${DEFAULT-VALUE}).")
    private int sequenceContextLength = 5;

    public static void main(String[] args) {
        int exitCode = new CommandLine(new AddAssortedAnnotationsToVcf()).execute(args);
        System.exit(exitCode);
    }

    @Override
    public Integer call() throws Exception {
        logger.info(getClass().getName() + " (" + getPackageNameAndVersion() + ")");

        IOUtil.assertFileIsReadable(inputVcfFile);
        IOUtil.assertFileIsReadable(referenceSequenceFile);
        IOUtil.assertFileIsWritable(outputVcfFile);

        // reference sequence
        ReferenceSequenceFile referenceSequence = ReferenceSequenceFileFactory
                .getReferenceSequenceFile(referenceSequenceFile);
        if (!referenceSequence.isIndexed()) {
            logger.error("Reference sequence FASTA file is not indexed");
            return 1;
        }

        VCFFileReader reader = new VCFFileReader(inputVcfFile, false);
        VCFHeader header = reader.getFileHeader();

        VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputVcfFile)
                .setOutputFileType(OutputType.VCF).setReferenceDictionary(header.getSequenceDictionary())
                .clearOptions();
        VariantContextWriter writer = builder.build();

        addInfoHeaderLines(header);
        writer.writeHeader(header);

        for (VariantContext variant : reader) {
            addMultiallelicFlag(variant);
            addSequenceContext(variant, referenceSequence);
            addIndelLength(variant);
            addHomopolymerLength(variant, referenceSequence);
            writer.add(variant);
        }

        writer.close();
        reader.close();
        referenceSequence.close();

        return 0;
    }

    /**
     * Add header lines for the added INFO fields.
     *
     * @param header the VCF header
     */
    private void addInfoHeaderLines(VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(MULTIALLELIC, 1, VCFHeaderLineType.Flag,
                "Indicates the variant is multiallelic."));
        header.addMetaDataLine(new VCFInfoHeaderLine(FIVE_PRIME_SEQUENCE_CONTEXT, 1, VCFHeaderLineType.String,
                "The 5' prime sequence context."));
        header.addMetaDataLine(new VCFInfoHeaderLine(THREE_PRIME_SEQUENCE_CONTEXT, 1, VCFHeaderLineType.String,
                "The 3' prime sequence context."));
        header.addMetaDataLine(
                new VCFInfoHeaderLine(INDEL_LENGTH, 1, VCFHeaderLineType.Integer, "The length of the indel."));
        header.addMetaDataLine(new VCFInfoHeaderLine(HOMOPOLYMER_LENGTH, 1, VCFHeaderLineType.Integer,
                "Length of the longest homopolymer run overlapping or immediately adjacent to the variant position."));
    }

    /**
     * Adds an INFO flag if the variant is multiallelic.
     *
     * @param variant
     */
    private void addMultiallelicFlag(VariantContext variant) {
        if (!variant.isBiallelic()) {
            variant.getCommonInfo().putAttribute(MULTIALLELIC, true);
        }
    }

    /**
     * Adds INFO entries for 5' and 3' sequence context.
     *
     * @param variant the variant
     */
    private void addSequenceContext(VariantContext variant, ReferenceSequenceFile referenceSequence) {
        variant.getCommonInfo().putAttribute(FIVE_PRIME_SEQUENCE_CONTEXT,
                referenceSequence.getSubsequenceAt(variant.getContig(), variant.getStart() - sequenceContextLength,
                        variant.getStart() - 1).getBaseString());
        variant.getCommonInfo().putAttribute(THREE_PRIME_SEQUENCE_CONTEXT, referenceSequence
                .getSubsequenceAt(variant.getContig(), variant.getEnd() + 1, variant.getEnd() + sequenceContextLength)
                .getBaseString());
    }

    /**
     * Adds an INFO entry for the indel length.
     *
     * @param variant the variant
     */
    private void addIndelLength(VariantContext variant) {
        if (variant.isIndel()) {
            int length = variant.getAlternateAllele(0).length() - variant.getReference().length();
            variant.getCommonInfo().putAttribute(INDEL_LENGTH, length);
        }
    }

    /**
     * Adds an INFO entry for the length of the longest homopolymer run
     * overlapping or immediately adjacent to the variant position. This is
     * useful for flagging indels that may be sequencing artifacts in
     * homopolymer regions.
     *
     * The method examines the reference sequence surrounding the variant,
     * extending outward from the variant position in both directions to find
     * the longest run of identical bases.
     *
     * @param variant           the variant
     * @param referenceSequence the reference sequence
     */
    private void addHomopolymerLength(VariantContext variant, ReferenceSequenceFile referenceSequence) {
        // use a window around the variant to search for homopolymer runs
        int windowSize = 20;
        int start = Math.max(1, variant.getStart() - windowSize);
        int end = variant.getEnd() + windowSize;

        // clamp to contig length
        int contigLength = referenceSequence.getSequenceDictionary()
                .getSequence(variant.getContig()).getSequenceLength();
        end = Math.min(end, contigLength);

        String bases = referenceSequence.getSubsequenceAt(variant.getContig(), start, end)
                .getBaseString().toUpperCase();

        // find the longest homopolymer run in the window
        int maxRun = 1;
        int currentRun = 1;
        for (int i = 1; i < bases.length(); i++) {
            if (bases.charAt(i) == bases.charAt(i - 1) && bases.charAt(i) != 'N') {
                currentRun++;
                if (currentRun > maxRun) {
                    maxRun = currentRun;
                }
            } else {
                currentRun = 1;
            }
        }

        variant.getCommonInfo().putAttribute(HOMOPOLYMER_LENGTH, maxRun);
    }
}
