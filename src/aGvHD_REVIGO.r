

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.
# --------------------------------------------------------------------------

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0032352","positive regulation of hormone metabolic process", 0.063,-2.767,-5.442, 1.079,70.6449,0.870,0.000),
c("GO:0071104","response to interleukin-9", 0.017, 6.126, 4.502, 0.602,70.6449,0.897,0.000),
c("GO:0072203","cell proliferation involved in metanephros development", 0.058,-6.361, 2.952, 1.041,70.6449,0.831,0.000),
c("GO:2000679","positive regulation of transcription regulatory region DNA binding", 0.127,-2.209,-6.684, 1.362,27.1711,0.961,0.019),
c("GO:0072593","reactive oxygen species metabolic process", 1.373, 0.910, 0.939, 2.378, 5.4134,0.977,0.026),
c("GO:0043900","regulation of multi-organism process", 1.379, 1.749, 2.365, 2.380, 6.0901,0.953,0.029),
c("GO:0009612","response to mechanical stimulus", 1.183, 3.038, 5.107, 2.314,12.4250,0.913,0.055),
c("GO:0050886","endocrine process", 0.456, 0.109, 6.211, 1.903,13.4136,0.942,0.081),
c("GO:2000615","regulation of histone H3-K9 acetylation", 0.075, 0.766,-4.524, 1.146,64.2226,0.883,0.091),
c("GO:0090399","replicative senescence", 0.069,-1.866, 6.070, 1.114,50.4606,0.902,0.099),
c("GO:0010566","regulation of ketone biosynthetic process", 0.069,-0.119,-7.372, 1.114,47.0966,0.904,0.099),
c("GO:1904996","positive regulation of leukocyte adhesion to vascular endothelial cell", 0.035,-5.438,-3.311, 0.845,44.1531,0.841,0.106),
c("GO:1905331","negative regulation of morphogenesis of an epithelium", 0.104,-5.567, 2.936, 1.279,41.5558,0.806,0.110),
c("GO:0048771","tissue remodeling", 0.837,-3.309, 3.289, 2.164, 8.6152,0.914,0.111),
c("GO:0034330","cell junction organization", 1.466,-1.273, 0.217, 2.407, 6.4223,0.941,0.116),
c("GO:1901216","positive regulation of neuron death", 0.398,-3.802,-4.603, 1.845,11.3943,0.872,0.126),
c("GO:0042107","cytokine metabolic process", 0.635, 2.294,-4.039, 2.045, 8.8306,0.967,0.126),
c("GO:0007568","aging", 1.633,-3.010, 6.221, 2.453, 6.7495,0.905,0.133),
c("GO:0046677","response to antibiotic", 0.283, 4.619, 5.792, 1.699, 7.6314,0.923,0.149),
c("GO:1902893","regulation of pri-miRNA transcription from RNA polymerase II promoter", 0.167, 2.894,-6.796, 1.477,17.6612,0.919,0.153),
c("GO:0043124","negative regulation of I-kappaB kinase/NF-kappaB signaling", 0.340, 6.063,-2.574, 1.778,23.5483,0.840,0.153),
c("GO:0001503","ossification", 2.147,-5.021, 1.528, 2.572, 4.5754,0.906,0.157),
c("GO:0019218","regulation of steroid metabolic process", 0.444,-0.879,-7.434, 1.892,11.4870,0.900,0.167),
c("GO:1903426","regulation of reactive oxygen species biosynthetic process", 0.398, 1.557,-6.911, 1.845,11.0383,0.908,0.173),
c("GO:0021700","developmental maturation", 1.420,-2.923, 7.335, 2.393, 6.7409,0.906,0.174),
c("GO:0043011","myeloid dendritic cell differentiation", 0.110,-1.535, 6.180, 1.301,39.2472,0.745,0.191),
c("GO:0070920","regulation of production of small RNA involved in gene silencing by RNA", 0.087, 5.983, 0.455, 1.204,33.6404,0.810,0.191),
c("GO:0001501","skeletal system development", 2.810,-5.136, 5.343, 2.688, 4.4431,0.860,0.194),
c("GO:0048799","animal organ maturation", 0.133,-3.826, 5.301, 1.380,32.1113,0.872,0.203),
c("GO:0016579","protein deubiquitination", 1.679, 2.180,-3.423, 2.465, 8.2787,0.950,0.203),
c("GO:2000146","negative regulation of cell motility", 1.310,-3.001,-3.099, 2.358, 7.3153,0.898,0.208),
c("GO:0001101","response to acid chemical", 1.858, 5.081, 5.865, 2.509, 5.4176,0.910,0.220),
c("GO:0002831","regulation of response to biotic stimulus", 0.802, 6.022,-1.769, 2.146,12.7980,0.854,0.229),
c("GO:0010817","regulation of hormone levels", 2.764,-0.071,-0.739, 2.681, 4.2472,0.940,0.233),
c("GO:2000810","regulation of bicellular tight junction assembly", 0.092,-0.927,-2.730, 1.230,32.1113,0.854,0.246),
c("GO:0009636","response to toxic substance", 1.264, 5.015, 6.174, 2.342, 4.3076,0.913,0.260),
c("GO:0051090","regulation of sequence-specific DNA binding transcription factor activity", 2.164, 1.877,-6.853, 2.575, 4.9985,0.900,0.276),
c("GO:0032332","positive regulation of chondrocyte differentiation", 0.110,-6.042, 2.506, 1.301,37.1815,0.759,0.288),
c("GO:1900543","negative regulation of purine nucleotide metabolic process", 0.381, 1.026,-6.216, 1.826,25.2303,0.868,0.289),
c("GO:0071773","cellular response to BMP stimulus", 0.923, 6.342, 4.197, 2.207, 9.8118,0.875,0.303),
c("GO:0033002","muscle cell proliferation", 0.969,-2.988, 1.932, 2.228,12.0247,0.951,0.319),
c("GO:0006775","fat-soluble vitamin metabolic process", 0.242,-3.606,-6.687, 1.633,16.4290,0.951,0.321),
c("GO:0007249","I-kappaB kinase/NF-kappaB signaling", 1.535, 7.294,-1.482, 2.427,11.2582,0.879,0.326),
c("GO:0002726","positive regulation of T cell cytokine production", 0.081, 0.006, 2.616, 1.176,37.1815,0.683,0.327),
c("GO:0051092","positive regulation of NF-kappaB transcription factor activity", 0.785, 0.167,-6.320, 2.137, 9.4825,0.866,0.341),
c("GO:0034105","positive regulation of tissue remodeling", 0.156,-6.265,-0.421, 1.447,19.0932,0.827,0.343),
c("GO:0043618","regulation of transcription from RNA polymerase II promoter in response to stress", 0.664, 3.902,-5.897, 2.064, 9.8118,0.844,0.346),
c("GO:0060231","mesenchymal to epithelial transition", 0.115,-4.467, 5.830, 1.322,35.3224,0.856,0.361),
c("GO:0010862","positive regulation of pathway-restricted SMAD protein phosphorylation", 0.277, 4.631,-3.025, 1.690,23.5483,0.778,0.374),
c("GO:0070757","interleukin-35-mediated signaling pathway", 0.121, 6.892, 2.985, 1.322,64.2226,0.827,0.378),
c("GO:0070997","neuron death", 1.760, 2.035,-1.654, 2.486, 5.1192,0.942,0.392),
c("GO:0070672","response to interleukin-15", 0.035, 6.146, 4.978, 0.845,50.4606,0.892,0.396),
c("GO:0043620","regulation of DNA-templated transcription in response to stress", 0.692, 4.117,-5.752, 2.083, 9.3776,0.845,0.404),
c("GO:0060044","negative regulation of cardiac muscle cell proliferation", 0.058,-5.221, 2.521, 1.041,35.3224,0.744,0.406),
c("GO:0003158","endothelium development", 0.687,-5.021, 5.090, 2.079,11.3943,0.870,0.419),
c("GO:0070106","interleukin-27-mediated signaling pathway", 0.012, 6.486, 3.321, 0.477,64.2226,0.853,0.420),
c("GO:0002040","sprouting angiogenesis", 0.415,-5.741, 3.908, 1.863, 8.4605,0.848,0.428),
c("GO:2001240","negative regulation of extrinsic apoptotic signaling pathway in absence of ligand", 0.196, 6.083,-1.201, 1.544,22.0765,0.820,0.430),
c("GO:0035723","interleukin-15-mediated signaling pathway", 0.017, 6.724, 3.324, 0.602,54.3422,0.845,0.431),
c("GO:0071496","cellular response to external stimulus", 1.587, 3.580, 2.935, 2.441, 7.7754,0.916,0.434),
c("GO:0035666","TRIF-dependent toll-like receptor signaling pathway", 0.162, 4.955,-1.915, 1.462,36.5405,0.721,0.437),
c("GO:0031098","stress-activated protein kinase signaling cascade", 1.506, 6.550,-2.665, 2.418, 7.6788,0.853,0.440),
c("GO:0060537","muscle tissue development", 2.072,-4.642, 5.938, 2.556, 4.6234,0.868,0.446),
c("GO:0002697","regulation of immune effector process", 1.922, 1.387, 4.361, 2.524, 5.0702,0.780,0.450),
c("GO:0018108","peptidyl-tyrosine phosphorylation", 2.135, 1.829,-4.677, 2.569, 5.1341,0.934,0.451),
c("GO:0050673","epithelial cell proliferation", 2.031,-2.840, 2.183, 2.548, 4.2353,0.948,0.453),
c("GO:2000352","negative regulation of endothelial cell apoptotic process", 0.162,-2.469,-4.179, 1.462,21.4075,0.909,0.458),
c("GO:0051492","regulation of stress fiber assembly", 0.421,-1.012,-2.486, 1.869,12.7672,0.868,0.465),
c("GO:0070102","interleukin-6-mediated signaling pathway", 0.058, 6.591, 2.938, 1.041,24.3603,0.836,0.467),
c("GO:0061756","leukocyte adhesion to vascular endothelial cell", 0.138,-5.429,-2.564, 1.398,16.0557,0.890,0.474),
c("GO:0001952","regulation of cell-matrix adhesion", 0.554,-5.952,-2.724, 1.987,12.2861,0.860,0.479),
c("GO:0034340","response to type I interferon", 0.485, 5.760, 2.830, 1.929,14.8726,0.743,0.483),
c("GO:0071305","cellular response to vitamin D", 0.110, 5.447, 4.390, 1.301,30.7152,0.866,0.485),
c("GO:0001773","myeloid dendritic cell activation", 0.167, 0.985, 3.922, 1.477,25.2303,0.840,0.488),
c("GO:0031589","cell-substrate adhesion", 1.783,-4.556,-2.010, 2.491, 7.7268,0.906,0.489),
c("GO:0070498","interleukin-1-mediated signaling pathway", 0.115, 6.720, 2.838, 1.322,13.9891,0.826,0.491),
c("GO:2000027","regulation of organ morphogenesis", 1.454,-5.431, 3.861, 2.403, 7.0929,0.764,0.496),
c("GO:0031664","regulation of lipopolysaccharide-mediated signaling pathway", 0.115, 6.128, 1.768, 1.322,28.2580,0.803,0.497),
c("GO:0003159","morphogenesis of an endothelium", 0.115,-5.442, 4.617, 1.322,37.1815,0.857,0.499));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );

