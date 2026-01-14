<!--
NOTES FOR CUSTOMIZATION:
- Replace [PROJECT_ID] with actual project number
    `perl -pe 's/\[PROJECT_ID\]/PROJECT_NO/g' template.md`
  eg: PROJECT_NO==12345_B
- Replace [Recipient Name] with recipient's name
- Replace [X.X] in resolution field with actual resolution used
- Include ribosomal gene paragraph only when applicable
- Can add/remove bullet points based on specific project parameters
-->

**Subject:** scRNA-seq Analysis Results

The results for your scRNA-seq project [PROJECT_ID] are now available. You can access them at:

  https://bicdelivery.mskcc.org/project/[PROJECT_ID]/seurat/r_001

Please find attached a detailed description of the output files and analyses (results.pdf).

Analysis Summary:

- Pipeline: Seurat v4 (current standard)
- QC filters: Standard cell filtering thresholds applied
- Clustering: Preliminary analysis performed at multiple resolutions (0.1, 0.2, 0.5, and 0.8)
- Initial differential expression analysis: Conducted at resolution **[X.X]** to identify cluster-specific markers
- Pathway analysis: Included for differentially expressed genes

<!-- INCLUDE ONLY IF RIBOSOMAL FILTERING WAS APPLIED: -->
For this analysis, we removed ribosomal genes as they appeared to be creating issues in the pathway analysis. If ribosomal genes are important to your experiment, please let us know and we can provide an alternative analysis.
<!-- END OPTIONAL SECTION -->

Next Steps:

Please review the clustering results and let us know if a different resolution would be more appropriate for your biological question. We can readily rerun the differential expression and pathway analyses at your preferred resolution.

If you have any questions or would like to discuss the results, please do not hesitate to reach out.

Best regards,

Nicholas D. Socci, PhD  
Director, Bioinformatics Core  
Memorial Sloan Kettering Cancer Center  
New York, NY 10021  
646.888.2608  
soccin@mskcc.org