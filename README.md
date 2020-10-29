# preciseTADworkshop

Precise prediction of TAD/loop boundaries, https://doi.org/10.1101/2020.09.03.282186


# European Bioconductor Conference 2020
https://openreview.net/group?id=bioconductor.org/EuroBioC/2020/Conference
Abstract submission for talks, posters, and workshops.

Virtual  Dec 14 2020  https://eurobioc2020.bioconductor.org/  eurobioc2020@stat.unipd.it

Please see the venue website for more information.

Submission Start: Aug 31 2020 12:00AM UTC-0, End: Oct 30 2020 12:00AM UTC-0

Add: 

**Title**

<!--A workshop for leveraging functional genomic annotations to--> preciseTAD: a machine-learning framework for predicting boundaries of 3D genomic elements

**Authors**

Search author profile by first, middle and last name or email address. If the profile is not found, you can add the author completing first, middle, last and name and author email address.

Spiro Stilianoudakis stilianoudasc@mymail.vcu.edu  
Mikhail G. Dozmorov mikhail.dozmorov@vcuhealth.org

**Keywords**

Hi-C, TAD, machine learning, classification, random forest, clustering, DBSCAN

**TL;DR**

preciseTAD, a supervised machine learning framework for precise prediction of boundaries of topologically associating domains (TADs) and chromatin loops. Our pre-trained models predict base-level domain boundaries using only genome annotation data.

**Abstract**

High-throughput chromosome conformation capture technology (Hi-C) revealed extensive DNA folding into discrete 3D structures referred to as Topologically Associating Domains (TADs) and loops.  TADs are critical for cellular processes like gene regulation and cell differentiation. The relatively low resolution of Hi-C data (tens of kilobases in size) prevents precise mapping of TAD boundaries by conventional TAD-callers. In contrast, the enrichment of high-resolution genomic annotations at boundaries (ChIP-seq, ~100 bases) offers a computational approach toward improved boundary identification. The workshop describes how to transform TAD-calling into a supervised machine learning framework. We demonstrate how to build and optimize a random forest model that prioritizes known molecular signatures of domain boundaries to predict each base's probability of being a boundary. We utilize density-based clustering (DBSCAN) and partitioning around medoids (PAM) to identify boundary regions and the most probable points. Boundaries identified by preciseTAD show strong signal enrichment of known boundary signatures compared with Arrowhead- and Peakachu boundaries. This workshop uses the R programming environment, the preciseTAD R package (https://bioconductor.org/packages/preciseTAD/), and can be performed with most operating systems on a single computer.

**Affiliation**

List your affiliation

**Gender**

**Bioconductor Package Maintenance**

preciseTAD

**Associated Packages**

preciseTAD

**Submission Type**

- Short workshop (30 - 45 minutes) such as a demonstration of a package or a workflow

**Submission Is Not Registration**

I am aware that this form does not constitute a registration. See the conference page for details about the registration.

**Code Of Conduct**

I agree to abide by the Bioconductor Code of Conduct (https://eurobioc2020.bioconductor.org/code_of_conduct) during the conference and in all interactions with the Bioconductor community. This includes prohibitions on harassment, discrimination, and intellectual bullying. It also includes a social media policy that anything presented at the conference may be posted publicly unless explicitly requested otherwise by the presenter.
