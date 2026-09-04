// =======================================================
// CytoOrg
// Quantifying cytoskeletal organization through skeleton-based image analysis
// Fiji Groovy script
//
// Version: 1.0
// Author: Perina Edoardo
// Date: 08/2026
//
// Run this single file from Fiji's Script Editor as Groovy.
// =======================================================



import ij.Executer
import ij.IJ
import ij.ImagePlus
import ij.Macro
import ij.WindowManager
import ij.CommandListener
import ij.measure.Calibration
import ij.measure.ResultsTable
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_
import sc.fiji.analyzeSkeleton.Edge
import sc.fiji.analyzeSkeleton.Graph
import sc.fiji.analyzeSkeleton.Point
import sc.fiji.analyzeSkeleton.SkeletonResult
import sc.fiji.analyzeSkeleton.Vertex

final class CytoOrgFilamentCandidate {
    int skeletonId
    double length
    double averageIntensity = Double.NaN
    Point startPoint
    Point endPoint
    boolean endpointsDefined

    CytoOrgFilamentCandidate(int skeletonId, double length,
                             Point startPoint, Point endPoint,
                             boolean endpointsDefined) {
        this.skeletonId = skeletonId
        this.length = length
        this.startPoint = startPoint
        this.endPoint = endPoint
        this.endpointsDefined = endpointsDefined
    }
}

final class CytoOrgMainPathResolution {
    double length
    Point startPoint
    Point endPoint
    boolean endpointsDefined
    Set<Edge> edges

    CytoOrgMainPathResolution(double length,
                              Point startPoint,
                              Point endPoint,
                              boolean endpointsDefined,
                              Set<Edge> edges) {
        this.length = length
        this.startPoint = startPoint
        this.endPoint = endPoint
        this.endpointsDefined = endpointsDefined
        this.edges = edges
    }
}

final class CytoOrgFilamentBridge implements CommandListener {
    static final String COMMAND = "CytoOrg Filament Bridge"
    static final String FILAMENT_TABLE = "CytoOrg Filaments v1.0"
    static final String ORIENTATION_TABLE = "CytoOrg Orientation v1.0"

    @Override
    String commandExecuting(String command) {
        if (command != COMMAND) {
            return command
        }

        String arguments = Macro.getOptions()

        try {
            createFilamentTables(arguments == null ? "" : arguments)
        }
        catch (Throwable error) {
            String detail = error.getMessage()

            if (detail == null || detail.trim().isEmpty()) {
                detail = error.getClass().getSimpleName()
            }

            IJ.error("CytoOrg Analyzer v1.0",
                    "Filament analysis failed:\n" + detail)
        }

        return null
    }

    private static void createFilamentTables(String arguments) {
        boolean useMainPathDefinition = arguments.contains("definition=main")
        boolean excludeSideBranches = arguments.contains("exclude_branches=true")
        boolean filterIntensity = arguments.contains("filter_intensity=true")
        ImagePlus sourceImage = WindowManager.getCurrentImage()

        if (sourceImage == null) {
            throw new IllegalStateException("No active skeleton image.")
        }

        ImagePlus analysisCopy = sourceImage.duplicate()
        analysisCopy.setTitle("CytoOrg_v1.0_filament_analysis_copy")

        try {
            AnalyzeSkeleton_ analyzer = new AnalyzeSkeleton_()
            analyzer.setup("", analysisCopy)

            SkeletonResult skeletonResult = analyzer.run(
                    AnalyzeSkeleton_.NONE,
                    false,
                    true,
                    analysisCopy,
                    true,
                    false)

            Graph[] graphs = skeletonResult.getGraph()
            List pathLengths = skeletonResult.getShortestPathList()
            def paths = analyzer.getShortestPathPoints()

            if (graphs == null || pathLengths == null || paths == null ||
                    graphs.length != pathLengths.size() ||
                    graphs.length != paths.length) {
                throw new IllegalStateException(
                        "AnalyzeSkeleton returned inconsistent graph and longest-path data.")
            }

            List<CytoOrgFilamentCandidate> selectedFilaments = []
            List<CytoOrgFilamentCandidate> orientationFilaments = []

            for (int skeletonIndex = 0;
                 skeletonIndex < graphs.length;
                 skeletonIndex++) {
                List path = paths[skeletonIndex]
                double reportedMainPathLength =
                        ((Number) pathLengths.get(skeletonIndex)).doubleValue()
                CytoOrgMainPathResolution mainPath = resolveMainPath(
                        graphs[skeletonIndex], reportedMainPathLength, path)
				
				boolean mainIntensityValid = mainPath.edges.every { Edge e ->
 					e.getColor() != 0
				}

				CytoOrgFilamentCandidate mainFilament = null

				if (mainPath.endpointsDefined && mainPath.length > 0 && (!filterIntensity || mainIntensityValid)) {
					mainFilament = new CytoOrgFilamentCandidate(
            		skeletonIndex + 1,
            		mainPath.length,
            		mainPath.startPoint,
            		mainPath.endPoint,
            		true)

    				mainFilament.averageIntensity = mainIntensityValid ? 255.0 : 0.0
				}
				
                List<CytoOrgFilamentCandidate> allSegments = []
                List<CytoOrgFilamentCandidate> mainPathSegments = []
                List<CytoOrgFilamentCandidate> sideSegments = []
                List edges = graphs[skeletonIndex].getEdges()

                for (Object edgeObject : edges) {
                    Edge edge = (Edge) edgeObject
                    if (filterIntensity && edge.getColor() == 0) {
    					continue
					}
                    Point segmentStartPoint = firstVertexPoint(edge.getV1())
                    Point segmentEndPoint = firstVertexPoint(edge.getV2())
                    boolean endpointsDefined =
                            segmentStartPoint != null && segmentEndPoint != null
                    CytoOrgFilamentCandidate segment =
                            new CytoOrgFilamentCandidate(
                                    skeletonIndex + 1,
                                    edge.getLength(),
                                    segmentStartPoint,
                                    segmentEndPoint,
                                    endpointsDefined)
					
					segment.averageIntensity = edge.getColor()
					
                    allSegments.add(segment)

                    if (mainPath.edges.contains(edge)) {
                        mainPathSegments.add(segment)
                    }
                    else {
                        sideSegments.add(segment)
                    }
                }

                if (useMainPathDefinition) {
                    if (mainFilament != null) {
                        selectedFilaments.add(mainFilament)
                    }

                    selectedFilaments.addAll(sideSegments)

                    if (mainFilament != null) {
                        orientationFilaments.add(mainFilament)
                    }

                    if (!excludeSideBranches) {
                        orientationFilaments.addAll(sideSegments)
                    }
                }
                else {
                    selectedFilaments.addAll(allSegments)

                    if (excludeSideBranches) {
                        orientationFilaments.addAll(mainPathSegments)
                    }
                    else {
                        orientationFilaments.addAll(allSegments)
                    }
                }
            }

            ResultsTable filamentTable = new ResultsTable()
            ResultsTable orientationTable = new ResultsTable()
            Calibration calibration = sourceImage.getCalibration()

            for (CytoOrgFilamentCandidate candidate : selectedFilaments) {
                addCandidateRow(filamentTable, candidate, calibration)
            }

            for (CytoOrgFilamentCandidate candidate : orientationFilaments) {
                addCandidateRow(orientationTable, candidate, calibration)
            }

            filamentTable.show(FILAMENT_TABLE)
            orientationTable.show(ORIENTATION_TABLE)
        }
        finally {
            analysisCopy.changes = false
            analysisCopy.close()
        }
    }

    private static CytoOrgMainPathResolution resolveMainPath(
            Graph graph,
            double reportedLength,
            List path) {
        Set<Edge> reportedEdges = newIdentityEdgeSet()

        if (path != null && path.size() >= 2 &&
                reportedLength > 0 && Double.isFinite(reportedLength)) {
            Set<String> pathPointKeys = new HashSet<String>()
            Set<String> pathLinkKeys = new HashSet<String>()

            for (int pointIndex = 0;
                 pointIndex < path.size();
                 pointIndex++) {
                Point pathPoint = (Point) path.get(pointIndex)
                pathPointKeys.add(pointKey(pathPoint))

                if (pointIndex > 0) {
                    pathLinkKeys.add(linkKey(
                            (Point) path.get(pointIndex - 1), pathPoint))
                }
            }

            for (Object edgeObject : graph.getEdges()) {
                Edge edge = (Edge) edgeObject

                if (edgeBelongsToMainPath(
                        edge, pathPointKeys, pathLinkKeys)) {
                    reportedEdges.add(edge)
                }
            }

            Point startPoint = (Point) path.get(0)
            Point endPoint = (Point) path.get(path.size() - 1)
            double reconstructedLength = 0

            for (Edge edge : reportedEdges) {
                reconstructedLength += edge.getLength()
            }

            boolean distinctEndpoints =
                    pointKey(startPoint) != pointKey(endPoint)
            double tolerance = 0.000000001 *
                    Math.max(1.0, Math.abs(reportedLength))

            if (distinctEndpoints && !reportedEdges.isEmpty() &&
                    Math.abs(reconstructedLength - reportedLength) <= tolerance) {
                return new CytoOrgMainPathResolution(
                        reportedLength,
                        startPoint,
                        endPoint,
                        true,
                        reportedEdges)
            }
        }

        return reconstructMainPathBetweenDistinctVertices(graph)
    }

    private static CytoOrgMainPathResolution
            reconstructMainPathBetweenDistinctVertices(Graph graph) {
        List vertices = graph.getVertices()
        List edges = graph.getEdges()
        int vertexCount = vertices.size()
        Set<Edge> emptyEdges = newIdentityEdgeSet()

        if (vertexCount < 2) {
            return new CytoOrgMainPathResolution(
                    0, null, null, false, emptyEdges)
        }

        double[][] distances = new double[vertexCount][vertexCount]
        int[][] nextVertex = new int[vertexCount][vertexCount]
        Edge[][] directEdges = new Edge[vertexCount][vertexCount]

        for (int firstIndex = 0;
             firstIndex < vertexCount;
             firstIndex++) {
            for (int secondIndex = 0;
                 secondIndex < vertexCount;
                 secondIndex++) {
                distances[firstIndex][secondIndex] =
                        firstIndex == secondIndex ?
                                0 : Double.POSITIVE_INFINITY
                nextVertex[firstIndex][secondIndex] = -1
            }
        }

        for (Object edgeObject : edges) {
            Edge edge = (Edge) edgeObject
            int firstIndex = vertices.indexOf(edge.getV1())
            int secondIndex = vertices.indexOf(edge.getV2())
            double edgeLength = edge.getLength()

            if (firstIndex < 0 || secondIndex < 0 ||
                    firstIndex == secondIndex || edgeLength <= 0 ||
                    !Double.isFinite(edgeLength)) {
                continue
            }

            if (edgeLength < distances[firstIndex][secondIndex]) {
                distances[firstIndex][secondIndex] = edgeLength
                distances[secondIndex][firstIndex] = edgeLength
                nextVertex[firstIndex][secondIndex] = secondIndex
                nextVertex[secondIndex][firstIndex] = firstIndex
                directEdges[firstIndex][secondIndex] = edge
                directEdges[secondIndex][firstIndex] = edge
            }
        }

        for (int intermediateIndex = 0;
             intermediateIndex < vertexCount;
             intermediateIndex++) {
            for (int firstIndex = 0;
                 firstIndex < vertexCount;
                 firstIndex++) {
                if (!Double.isFinite(
                        distances[firstIndex][intermediateIndex])) {
                    continue
                }

                for (int secondIndex = 0;
                     secondIndex < vertexCount;
                     secondIndex++) {
                    if (!Double.isFinite(
                            distances[intermediateIndex][secondIndex])) {
                        continue
                    }

                    double candidateDistance =
                            distances[firstIndex][intermediateIndex] +
                            distances[intermediateIndex][secondIndex]

                    if (candidateDistance <
                            distances[firstIndex][secondIndex]) {
                        distances[firstIndex][secondIndex] = candidateDistance
                        nextVertex[firstIndex][secondIndex] =
                                nextVertex[firstIndex][intermediateIndex]
                    }
                }
            }
        }

        double longestShortestDistance = 0
        int pathStartIndex = -1
        int pathEndIndex = -1

        for (int firstIndex = 0;
             firstIndex < vertexCount;
             firstIndex++) {
            for (int secondIndex = firstIndex + 1;
                 secondIndex < vertexCount;
                 secondIndex++) {
                double distance = distances[firstIndex][secondIndex]

                if (Double.isFinite(distance) &&
                        distance > longestShortestDistance) {
                    longestShortestDistance = distance
                    pathStartIndex = firstIndex
                    pathEndIndex = secondIndex
                }
            }
        }

        if (pathStartIndex < 0 || pathEndIndex < 0) {
            return new CytoOrgMainPathResolution(
                    0, null, null, false, emptyEdges)
        }

        Set<Edge> reconstructedEdges = newIdentityEdgeSet()
        int currentIndex = pathStartIndex
        int reconstructionSteps = 0

        while (currentIndex != pathEndIndex) {
            int followingIndex = nextVertex[currentIndex][pathEndIndex]

            if (followingIndex < 0 ||
                    reconstructionSteps >= vertexCount) {
                return new CytoOrgMainPathResolution(
                        0, null, null, false, emptyEdges)
            }

            Edge pathEdge = directEdges[currentIndex][followingIndex]

            if (pathEdge == null) {
                return new CytoOrgMainPathResolution(
                        0, null, null, false, emptyEdges)
            }

            reconstructedEdges.add(pathEdge)
            currentIndex = followingIndex
            reconstructionSteps++
        }

        Point startPoint = firstVertexPoint(
                (Vertex) vertices.get(pathStartIndex))
        Point endPoint = firstVertexPoint(
                (Vertex) vertices.get(pathEndIndex))
        boolean endpointsDefined = startPoint != null && endPoint != null &&
                pointKey(startPoint) != pointKey(endPoint)

        if (!endpointsDefined || reconstructedEdges.isEmpty()) {
            return new CytoOrgMainPathResolution(
                    0, null, null, false, emptyEdges)
        }

        return new CytoOrgMainPathResolution(
                longestShortestDistance,
                startPoint,
                endPoint,
                true,
                reconstructedEdges)
    }

    private static Set<Edge> newIdentityEdgeSet() {
        return Collections.newSetFromMap(
                new IdentityHashMap<Edge, Boolean>())
    }

    private static String pointKey(Point point) {
        return point.x + ":" + point.y + ":" + point.z
    }

    private static String linkKey(Point firstPoint, Point secondPoint) {
        String firstKey = pointKey(firstPoint)
        String secondKey = pointKey(secondPoint)

        if (firstKey.compareTo(secondKey) <= 0) {
            return firstKey + "|" + secondKey
        }

        return secondKey + "|" + firstKey
    }

    private static boolean edgeBelongsToMainPath(
            Edge edge,
            Set<String> pathPointKeys,
            Set<String> pathLinkKeys) {
        List slabs = edge.getSlabs()

        if (slabs != null && !slabs.isEmpty()) {
            for (Object slabObject : slabs) {
                if (!pathPointKeys.contains(pointKey((Point) slabObject))) {
                    return false
                }
            }

            return true
        }

        List firstVertexPoints = edge.getV1().getPoints()
        List secondVertexPoints = edge.getV2().getPoints()

        for (Object firstObject : firstVertexPoints) {
            for (Object secondObject : secondVertexPoints) {
                if (pathLinkKeys.contains(linkKey(
                        (Point) firstObject, (Point) secondObject))) {
                    return true
                }
            }
        }

        return false
    }

    private static Point firstVertexPoint(Vertex vertex) {
        List points = vertex.getPoints()
        return points == null || points.isEmpty() ?
                null : (Point) points.get(0)
    }

    private static void addCandidateRow(
            ResultsTable table,
            CytoOrgFilamentCandidate candidate,
            Calibration calibration) {
        double startX = Double.NaN
        double startY = Double.NaN
        double startZ = Double.NaN
        double endX = Double.NaN
        double endY = Double.NaN
        double endZ = Double.NaN
        double euclideanDistance = Double.NaN

        if (candidate.endpointsDefined) {
            startX = candidate.startPoint.x * calibration.pixelWidth
            startY = candidate.startPoint.y * calibration.pixelHeight
            startZ = candidate.startPoint.z * calibration.pixelDepth
            endX = candidate.endPoint.x * calibration.pixelWidth
            endY = candidate.endPoint.y * calibration.pixelHeight
            endZ = candidate.endPoint.z * calibration.pixelDepth

            double deltaX = endX - startX
            double deltaY = endY - startY
            double deltaZ = endZ - startZ
            euclideanDistance = Math.sqrt(
                    deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ)
        }

        table.incrementCounter()
        table.addValue("Skeleton ID", candidate.skeletonId)
        table.addValue("Filament length", candidate.length)
        table.addValue("Euclidean distance", euclideanDistance)
        table.addValue("V1 x", startX)
        table.addValue("V1 y", startY)
        table.addValue("V1 z", startZ)
        table.addValue("V2 x", endX)
        table.addValue("V2 y", endY)
        table.addValue("V2 z", endZ)
        table.addValue("average intensity", candidate.averageIntensity)
    }
}

CytoOrgFilamentBridge bridge = new CytoOrgFilamentBridge()
Executer.addCommandListener(bridge)

try {
    String macroSource = '''

dir = getDirectory("Choose input folder");

micro = fromCharCode(181);
sup2 = fromCharCode(178);
helpUrl = "https://github.com/EP-bioimage-tools/CytoOrg-ImageJ";

filamentDefinitions = newArray(
    "Segments between junctions/endpoints",
    "Main filament (longest shortest path)"
);

loopPruningMethods = newArray(
    "None",
    "Shortest branch"
);

Dialog.create("CytoOrg Analyzer");

Dialog.addMessage("Select the analyses and basic settings.\\nEnable Advanced options to modify the defaults after clicking OK.", 15, "#555555");

Dialog.setInsets(15, 0, 0);
Dialog.addMessage("GENERAL SETTINGS", 15, "#1F4E79");

Dialog.addCheckbox("Use segcell masks for multiple cells", false);
Dialog.addMessage("For every original image, add a mask named\\nsegcell_<complete image filename> to the input folder.", 12, "#777777");

Dialog.setInsets(10, 0, 4);
Dialog.addNumber("Pixel size", 0.0547606, 3, 9, " " + micro + "m/pixel");

Dialog.setInsets(10, 0, 4);
Dialog.addChoice("Filament definition", filamentDefinitions, filamentDefinitions[0]);
Dialog.addMessage("Segments: each section between junctions or end-points is one filament.\\nMain filament: the complete longest shortest path is one filament; off-path\\nsegments remain side branches. This affects count, length and orientation.", 12, "#777777");

Dialog.setInsets(20, 0, 0);
Dialog.addMessage("FILAMENT THICKNESS", 15, "#1F4E79");

Dialog.addCheckbox("Measure filament thickness", true);
Dialog.setInsets(8, 0, 4);
Dialog.addNumber("Minimum thickness", 0.2, 2, 5, " " + micro + "m");
Dialog.addCheckbox("Advanced thickness settings", false);

Dialog.setInsets(25, 0, 0);
Dialog.addMessage("SKELETON ANALYSIS", 15, "#1F4E79");

Dialog.setInsets(8, 0, 4);
Dialog.addNumber("Tubeness sigma", 0.3, 2, 5, " " + micro + "m");

Dialog.addCheckbox("Estimate sigma from filament thickness", false);
Dialog.addMessage("Thickness measurement is required for automatic sigma.", 12, "#777777");

Dialog.setInsets(8, 0, 4);
Dialog.addNumber("Minimum length", 0.2, 2, 5, " " + micro + "m");

Dialog.addCheckbox("Advanced skeleton settings", false);

Dialog.setInsets(20, 0, 0);
Dialog.addMessage("FILAMENT ORIENTATION ANALYSIS", 15, "#1F4E79");

Dialog.addCheckbox("Analyze filament orientation", true);
Dialog.addCheckbox("Advanced orientation settings", false);

Dialog.setInsets(25, 0, 0);
Dialog.addMessage("OUTPUT", 15, "#1F4E79");

Dialog.addString("Output folder name", "Results", 20);

Dialog.addHelp(helpUrl);
Dialog.show();

MultiCellAnalysis = Dialog.getCheckbox();

pixelSize = Dialog.getNumber();
filamentDefinition = Dialog.getChoice();

measureThickness = Dialog.getCheckbox();
minBranchThickness = Dialog.getNumber();
ThicknessOptions = Dialog.getCheckbox();

sigmaTubeness = Dialog.getNumber();
automaticSigma = Dialog.getCheckbox();
minBranchLength = Dialog.getNumber();
SkeletonOptions = Dialog.getCheckbox();

OrientationAnalysis = Dialog.getCheckbox();
OrientationOptions = Dialog.getCheckbox();

OutputFolderName = Dialog.getString();

useMainPathDefinition = false;

if (filamentDefinition == filamentDefinitions[1]) {
    useMainPathDefinition = true;
}

if (automaticSigma && !measureThickness) {
    exit("Automatic sigma estimation requires filament thickness measurement.");
}

thicknessMaxSuggested = minBranchThickness * 7;

if (measureThickness && ThicknessOptions) {

    Dialog.create("Filament Thickness");

    Dialog.addMessage("DETECTION", 15, "#1F4E79");

    Dialog.setInsets(8, 0, 4);
    Dialog.addSlider("Intensity threshold (0 = automatic)", 0, 255, 40);
 
    Dialog.setInsets(18, 0, 4);
    Dialog.addMessage("THICKNESS MAP DISPLAY", 15, "#1F4E79");
    Dialog.addMessage("Controls only the color scale of saved thickness maps.\\nMeasured values are not clipped or excluded.", 12, "#777777");

    Dialog.setInsets(8, 0, 4);
    Dialog.addNumber("Color scale minimum", 0, 2, 5, " " + micro + "m");

    Dialog.setInsets(8, 0, 4);
    Dialog.addNumber("Color scale maximum", thicknessMaxSuggested, 2, 5, " " + micro + "m");

    Dialog.addHelp(helpUrl);
    Dialog.show();

    IntensityThreshold = round(Dialog.getNumber());
    thicknessMin = Dialog.getNumber();
    thicknessMax = Dialog.getNumber();

} else {

    IntensityThreshold = 0;
    thicknessMin = 0;
    thicknessMax = thicknessMaxSuggested;
}

ShortIntermediateSuggested = minBranchLength * 9;
IntermediateLongSuggested = (minBranchLength * 18) - minBranchLength;

if (SkeletonOptions) {

    Dialog.create("Skeleton Analysis");

    Dialog.setInsets(12, 0, 2);
    Dialog.addMessage("PRUNING", 15, "#1F4E79");

    Dialog.setInsets(8, 0, 4);
    Dialog.addChoice("Loop pruning method", loopPruningMethods, loopPruningMethods[0]);
    Dialog.addMessage("Shortest branch cuts the shortest branch found in each detected loop.", 12, "#777777");

    Dialog.addCheckbox("Prune branches ending at end-points", false);
    Dialog.addMessage("Removes complete terminal branches.", 12, "#777777");

    Dialog.setInsets(18, 0, 2);
    Dialog.addMessage("SKELETON ANALYSIS", 15, "#1F4E79");

	Dialog.addCheckbox("Group filaments by length", false);
	Dialog.addMessage("Adds filament count and density results for the short, intermediate\\nand long groups to the main results table.", 12, "#777777");

    Dialog.setInsets(8, 0, 4);
    Dialog.addNumber("Short / intermediate threshold", ShortIntermediateSuggested, 2, 5, " " + micro + "m");

    Dialog.setInsets(8, 0, 4);
    Dialog.addNumber("Intermediate / long threshold", IntermediateLongSuggested, 2, 5, " " + micro + "m");

    Dialog.setInsets(12, 20, 0);
    Dialog.addCheckbox("Export individual filament data", false);

    Dialog.addCheckbox("Export grouped filament data", false);

    Dialog.addCheckbox("Save enhanced skeleton", true);
    Dialog.addMessage("Enhanced maps use thicker lines for easier visualization; measurements\\ncontinue to use the original skeleton.", 12, "#777777");

    Dialog.setInsets(18, 0, 2);
    Dialog.addMessage("BRANCHING", 15, "#1F4E79");

    Dialog.addCheckbox("Analyze filament branching", true);
    Dialog.addCheckbox("Save branching map", false);

    Dialog.addCheckbox("Save enhanced branching map", true);

    Dialog.addHelp(helpUrl);
    Dialog.show();

    loopPruningMethod = Dialog.getChoice();
    pruneEndBranches = Dialog.getCheckbox();

	StratifyLengths = Dialog.getCheckbox();
    ShortIntermediate = Dialog.getNumber();
    IntermediateLong = Dialog.getNumber();

    saveIndividualFilaments = Dialog.getCheckbox();
    saveGroupsFilaments = Dialog.getCheckbox();
    saveEnhancedSkeleton = Dialog.getCheckbox();

    AnalyzeBranching = Dialog.getCheckbox();
    saveBranchMap = Dialog.getCheckbox();
    saveEnhancedBranchMap = Dialog.getCheckbox();

    if (!AnalyzeBranching) {
        saveBranchMap = false;
        saveEnhancedBranchMap = false;
    }

} else {

    loopPruningMethod = loopPruningMethods[0];
    pruneEndBranches = false;

	StratifyLengths = false;
    ShortIntermediate = ShortIntermediateSuggested;
    IntermediateLong = IntermediateLongSuggested;

    saveIndividualFilaments = false;
    saveGroupsFilaments = false;
    saveEnhancedSkeleton = true;

    AnalyzeBranching = true;
    saveBranchMap = false;
    saveEnhancedBranchMap = true;
}

methods = newArray(
    "Cubic Spline",
    "Finite Difference",
    "Fourier",
    "Riesz Filters",
    "Gaussian",
    "Hessian"
);

if (OrientationAnalysis && OrientationOptions) {

    Dialog.create("Filament Orientation");

    Dialog.addMessage("BRANCH ORIENTATION AND ANGULAR DISTRIBUTION", 15, "#1F4E79");

    Dialog.addCheckbox("Exclude side branches from orientation statistics", false);
    Dialog.addMessage("Main-filament definition: keeps only complete longest paths.\\nSegment definition: keeps only the segments belonging to those paths.\\nWhen this option is disabled, side branches are also included.", 12, "#777777");

    Dialog.setInsets(8, 0, 4);
    Dialog.addNumber("Maximum tortuosity included", 1.5, 2, 6, "");

    Dialog.addCheckbox("Weight orientation statistics by filament length", true);

    Dialog.setInsets(10, 20, 2);
    Dialog.addCheckbox("Save angular distribution plot", true);
    Dialog.addMessage("The following settings are used only when the plot is saved.", 12, "#777777");

    Dialog.setInsets(8, 0, 4);
    Dialog.addNumber("Angular bin width", 5, 0, 5, " degrees");

    Dialog.setInsets(8, 0, 4);
    Dialog.addSlider("Radial-axis maximum (0 = automatic)", 0, 100, 0);
    Dialog.addMessage("0: automatic scale calculated separately for each plot.\\n1-100%: fixed maximum shared by all plots; higher values saturate\\n(10% generally gives a useful common scale across images).", 12, "#777777");

    Dialog.setInsets(20, 0, 2);
    Dialog.addMessage("DIRECTION-CHANGE ANALYSIS", 15, "#1F4E79");
    Dialog.addMessage("Measures how often the local filament direction changes by at least\\nthe selected angle.", 12, "#777777");

    Dialog.setInsets(8, 0, 4);
    Dialog.addNumber("Orientation window sigma", 2);

    Dialog.setInsets(8, 0, 4);
    Dialog.addChoice("Gradient method", methods, "Cubic Spline");

    Dialog.setInsets(8, 0, 4);
    Dialog.addSlider("Minimum direction change (degrees)", 1, 90, 45);

    Dialog.addCheckbox("Save orientation map", false);
    Dialog.addCheckbox("Save enhanced orientation map", true);

    Dialog.addHelp(helpUrl);
    Dialog.show();

    excludeSideBranchesFromOrientation = Dialog.getCheckbox();
    maxOrientationTortuosity = Dialog.getNumber();
    weightOrientationByLength = Dialog.getCheckbox();
    saveAngularDistribution = Dialog.getCheckbox();
    orientationBinWidth = Dialog.getNumber();
    angularPlotMaximum = Dialog.getNumber();
    orientationSigma = Dialog.getNumber();
    gradientMethod = Dialog.getChoice();
    directionChangeThreshold = Dialog.getNumber();
    saveOrientationMap = Dialog.getCheckbox();
    saveEnhancedOrientationMap = Dialog.getCheckbox();

} else {

    orientationSigma = 2;
    gradientMethod = "Cubic Spline";
    directionChangeThreshold = 45;
    maxOrientationTortuosity = 1.5;
    excludeSideBranchesFromOrientation = false;
    weightOrientationByLength = true;
    orientationBinWidth = 5;
    angularPlotMaximum = 0;
    saveOrientationMap = false;
    saveEnhancedOrientationMap = OrientationAnalysis;
    saveAngularDistribution = OrientationAnalysis;
}

if (OrientationAnalysis && maxOrientationTortuosity < 1) {
    exit("Maximum tortuosity must be at least 1.");
}

if (OrientationAnalysis &&
    (orientationBinWidth != floor(orientationBinWidth) ||
     orientationBinWidth < 1 || orientationBinWidth > 90 ||
     180 % orientationBinWidth != 0)) {
    exit("Angular bin width must be an integer divisor of 180 (for example 1, 2, 3, 5 or 10 degrees).");
}

if (OrientationAnalysis && (angularPlotMaximum < 0 || angularPlotMaximum > 100)) {
    exit("Radial-axis maximum must be between 0 and 100 percent.");
}

files = getFileList(dir);
originalFiles = newArray(0);

for (i = 0; i < files.length; i++) {
    name = files[i];
    lowerName = toLowerCase(name);

    if (isSupportedImageFile(name) && !startsWith(lowerName, "segcell_")) {
        originalFiles = Array.concat(originalFiles, name);
    }
}

if (originalFiles.length == 0) {
    exit("No valid original image files found in selected folder.");
}

if (MultiCellAnalysis) {

    missingMasks = "";

    for (i = 0; i < originalFiles.length; i++) {
        maskName = "segcell_" + originalFiles[i];

        if (!File.exists(dir + maskName)) {
            missingMasks += "\\n- " + maskName;
        }
    }

    if (missingMasks != "") {
        exit("Multi-cell analysis cannot start. Missing mask file(s):" + missingMasks);
    }
}

if (OutputFolderName=="") {
    OutputFolderName = "Results";
}

outputDir = dir + OutputFolderName + File.separator;

if (!File.exists(outputDir)) {
    File.makeDirectory(outputDir);
}

saveSettings;
setBatchMode(true);
setOption("BlackBackground", true);

analysisDir = dir;
temporaryDir = "";
validFiles = originalFiles;
roiControlFiles = newArray(0);
segmentedCellAreas = newArray(0);
segmentedCellRoundness = newArray(0);
segmentedCellSolidity = newArray(0);

if (MultiCellAnalysis) {

    for (i = 0; i < originalFiles.length; i++) {

        imageName = originalFiles[i];
        originalPath = dir + imageName;
        maskPath = dir + "segcell_" + imageName;

        dimensionError = getDimensionMismatchMessage(originalPath, maskPath, imageName);

        if (dimensionError != "") {
            setBatchMode(false);
            exit(dimensionError);
        }
    }

    temporaryDir = dir + ".cytoorg_temp_" + d2s(getTime(), 0) + File.separator;
    File.makeDirectory(temporaryDir);

    validFiles = newArray(0);
    createdCellImages = 0;

    for (i = 0; i < originalFiles.length; i++) {

        imageName = originalFiles[i];
        originalPath = dir + imageName;
        maskPath = dir + "segcell_" + imageName;
        baseName = getBaseName(imageName);

        print("Finding cells in mask " + (i+1) + "/" + originalFiles.length + ": segcell_" + imageName);

        roiManager("Reset");
        run("Clear Results");
		
		open(maskPath);

		setVoxelSize(pixelSize, pixelSize, 1, "um");

		run("8-bit");
		setThreshold(255, 255);

		run("Set Measurements...", "area perimeter bounding fit shape feret's redirect=None decimal=3");

		run("Analyze Particles...", "size=0-Infinity show=Nothing clear add");

		nCells = roiManager("count");

        if (nCells == 0) {
            close();
            roiManager("Reset");
            deleteTemporaryFolder(temporaryDir);
            setBatchMode(false);
            exit("No cell ROI was found in mask:\\nsegcell_" + imageName);
        }

		for (cellIndex = 0; cellIndex < nCells; cellIndex++) {

    		segmentedCellAreas = Array.concat(
        		segmentedCellAreas,
        		getResult("Area", cellIndex)
    		);

    		segmentedCellRoundness = Array.concat(
        		segmentedCellRoundness,
        		getResult("Round", cellIndex)
    		);

    		segmentedCellSolidity = Array.concat(
        		segmentedCellSolidity,
        		getResult("Solidity", cellIndex)
    		);
		}

        close();

        roiControlFileName = createRoiControlImage(originalPath, baseName, temporaryDir, nCells);
        roiControlFiles = Array.concat(roiControlFiles, roiControlFileName);

        for (cellIndex = 0; cellIndex < nCells; cellIndex++) {

            open(originalPath);
            roiManager("Select", cellIndex);

            setBackgroundColor(0, 0, 0);
            run("Clear Outside");
            run("Select None");

            cellImageName = baseName + "_cell" + (cellIndex+1) + ".tif";
            cellImagePath = temporaryDir + cellImageName;

            saveAs("Tiff", cellImagePath);
            close();

            validFiles = Array.concat(validFiles, cellImageName);
            createdCellImages++;
        }

        roiManager("Reset");
    }

    analysisDir = temporaryDir;
    run("Clear Results");

	print("Number of cells to be analyzed: " + createdCellImages);
}

if (validFiles.length == 0) {
    if (MultiCellAnalysis) deleteTemporaryFolder(temporaryDir);
    setBatchMode(false);
    exit("No valid image files are available for analysis.");
}

OrganizeResults = validFiles.length >= 10;

thicknessOutputDir = outputDir;
skeletonOutputDir = outputDir;
improvedSkeletonOutputDir = outputDir;
individualFilamentsOutputDir = outputDir;
groupedFilamentsOutputDir = outputDir;
branchingOutputDir = outputDir;
improvedBranchingOutputDir = outputDir;
orientationOutputDir = outputDir;
improvedOrientationOutputDir = outputDir;
angularDistributionOutputDir = outputDir;
roiControlOutputDir = outputDir;

if (OrganizeResults) {

    skeletonOutputDir = getCategoryOutputDirectory(outputDir, "skeleton");

    if (measureThickness)
        thicknessOutputDir = getCategoryOutputDirectory(outputDir, "thickness_map");

    if (saveEnhancedSkeleton)
        improvedSkeletonOutputDir = getCategoryOutputDirectory(outputDir, "skeleton_improved");

    if (saveIndividualFilaments)
        individualFilamentsOutputDir = getCategoryOutputDirectory(outputDir, "single_filaments");

    if (saveGroupsFilaments)
        groupedFilamentsOutputDir = getCategoryOutputDirectory(outputDir, "grouped_filaments");

    if (saveBranchMap)
        branchingOutputDir = getCategoryOutputDirectory(outputDir, "branching_map");

    if (saveEnhancedBranchMap)
        improvedBranchingOutputDir = getCategoryOutputDirectory(outputDir, "branching_map_improved");

    if (saveOrientationMap)
        orientationOutputDir = getCategoryOutputDirectory(outputDir, "orientation_map");

    if (saveEnhancedOrientationMap)
        improvedOrientationOutputDir = getCategoryOutputDirectory(outputDir, "orientation_map_improved");

    if (saveAngularDistribution)
        angularDistributionOutputDir = getCategoryOutputDirectory(outputDir, "angular_distribution");

    if (MultiCellAnalysis)
		roiControlOutputDir = getCategoryOutputDirectory(outputDir, "ROI_overlays");
}

if (MultiCellAnalysis) {

    roiCopyError = copyRoiControlImages(roiControlFiles, temporaryDir, roiControlOutputDir);

    if (roiCopyError != "") {
        deleteTemporaryFolder(temporaryDir);
        setBatchMode(false);
        exit(roiCopyError);
    }
}

resultRow = 0;
Table.create("Cytoskeletal Organization - Results");

for (i = 0; i < validFiles.length; i++) {

    imageName = validFiles[i];
    imagePath = analysisDir + imageName;
    
	baseName = getBaseName(imageName);
    
    print("Analyzing image " + (i+1) + "/" + validFiles.length + ": " + imageName);

	open(imagePath);
	run("Set Scale...", "distance=1 known=" + pixelSize + " unit=um global");

	if (MultiCellAnalysis) {

    	CellArea = segmentedCellAreas[i];
    	Roundness = segmentedCellRoundness[i];
    	Solidity = segmentedCellSolidity[i];

	} else {

    run("8-bit");

    setThreshold(1, 255);
    run("Convert to Mask");
    run("Fill Holes");

    run("Set Measurements...","area mean centroid center perimeter bounding fit shape feret's area_fraction redirect=None decimal=3");

    run("Clear Results");
    run("Analyze Particles...", "clear");

    CellArea = getResult("Area", 0);
    Perimeter = getResult("Perim.", 0);
    MajorAxis = getResult("Major", 0);
    MinorAxis = getResult("Minor", 0);
    AspectRatio = getResult("AR", 0);
    Roundness = getResult("Round", 0);
    Solidity = getResult("Solidity", 0);
    Angle = getResult("Angle", 0);
	
	}

	close();
		
	if (measureThickness) {
		
		minBranchThickness_px = minBranchThickness / pixelSize;
		
		open(imagePath);
		
		run("8-bit");
		
		if(IntensityThreshold == 0){
			
			setAutoThreshold("Otsu dark");
			run("Convert to Mask");
			run("Local Thickness (complete process)", "threshold=128");

		} else {
			
			run("Local Thickness (complete process)", "threshold=" + IntensityThreshold);
		
		}
		
		width = getWidth();
		height = getHeight();

		thicknessValues = newArray(width * height);

		sumThickness = 0;
		nPixels = 0;

		for (x = 0; x < width; x++) {
  		  for (y = 0; y < height; y++) {

       	 	value = getPixel(x, y);

        	if (value > minBranchThickness_px) {
        		
        		thicknessValues[nPixels] = value;
            	sumThickness += value;
            	nPixels++;
        		
        		} else {
        			
        		setPixel(x, y, 0);
        			
        		}
    		}
		}

		thicknessValues = Array.trim(thicknessValues, nPixels);
		Array.sort(thicknessValues);

		if (nPixels % 2 == 1) {

    		median_px = thicknessValues[(nPixels-1)/2];

		} else {

    		median_px = (thicknessValues[nPixels/2 - 1] + thicknessValues[nPixels/2]) / 2;
		}

		median_um = median_px * pixelSize;
		
		q1_px = thicknessValues[floor(0.25 * (nPixels - 1))];
		q3_px = thicknessValues[floor(0.75 * (nPixels - 1))];
		q1_um = q1_px * pixelSize;
		q3_um = q3_px * pixelSize;

		meanThickness_px = sumThickness / nPixels;
		meanThickness_um = meanThickness_px * pixelSize;
		
		thicknessMin_px = thicknessMin / pixelSize;
		thicknessMax_px = thicknessMax / pixelSize;
		
		setMinAndMax(thicknessMin_px, thicknessMax_px);
		run("RGB Color");
		
		saveAs("Tiff", thicknessOutputDir + "Thickness_" + baseName + ".tif");
		
		close("*");
	}

   		open(imagePath);
   		
   		run("8-bit");
   		
   		if (automaticSigma) {
   			
   			sigma = meanThickness_um/2;
			run("Tubeness", "sigma=" + sigma + " use");
			setOption("ScaleConversions", true);
			
		} else {
    
    		run("Tubeness", "sigma=" + sigmaTubeness + " use");
			setOption("ScaleConversions", true);
    
		}
   		
   		run("8-bit");
   		setThreshold(1, 255);
        run("Convert to Mask");
       
       	run("Create Selection");		
		getStatistics(area, mean, min, max, std, histogram);
		CytoskeletonArea = area;
    
    	run("Skeletonize");
    	
    	run("Clear Results");
		roiManager("Reset");
    	
    	run("Analyze Particles...", "size=0-infinity clear add");
    	
    	nROI = roiManager("count");

		for (p = nROI-1; p >= 0; p--) {

  		  	perimeterROI = getResult("Perim.", p);
  		  	AreaROI = getResult("Area", p);

    		if (perimeterROI < 2.1*minBranchLength || AreaROI < 0.05*minBranchLength) {

        	roiManager("Select", p);

        	setForegroundColor(0,0,0);
       		run("Fill");
       		
    		}
		}
    	
		roiManager("Reset");
		run("Remove Overlay");
		run("Select None");

    	saveAs("Tiff", skeletonOutputDir + "Skeleton_" + baseName + ".tif");
    	SkeletonPath = skeletonOutputDir + "Skeleton_" + baseName + ".tif";
		
		close("*");
    
    	run("Clear Results");
    	
    	open(SkeletonPath);
    	skeletonAnalysisImageID = getImageID();

        analyzeSkeletonOptions = "prune=none";

        if (loopPruningMethod == loopPruningMethods[1])
            analyzeSkeletonOptions = "prune=[shortest branch]";

        if (pruneEndBranches)
            analyzeSkeletonOptions += " prune_0";

        analyzeSkeletonOptions += " calculate show";

	    	run("Analyze Skeleton (2D/3D)", analyzeSkeletonOptions);

        selectImage(skeletonAnalysisImageID);

        if (loopPruningMethod != loopPruningMethods[0] || pruneEndBranches) {
            fileDeleteStatus = File.delete(SkeletonPath);
            saveAs("Tiff", SkeletonPath);
        }

        if (saveEnhancedSkeleton) {
            run("Duplicate...", "title=CytoOrg_Enhanced_Skeleton_Temp");
            run("Maximum...", "radius=3");
            saveAs("Tiff", improvedSkeletonOutputDir + "Improved_Skeleton_" + baseName + ".tif");
            close();
            selectImage(skeletonAnalysisImageID);
        }
    
        totalJunctions = 0;
        analyzeSkeletonResultsSize = Table.size("Results");

        for (a = 0; a < analyzeSkeletonResultsSize; a++) {
            AverageBranchLength = Table.get("Average Branch Length", a, "Results");

            if (AverageBranchLength > minBranchLength) {
                totalJunctions += Table.get("# Junctions", a, "Results");
            }
        }

        nSkeletonSegments = 0;
        branchInformationSize = Table.size("Branch information");

        for (b = 0; b < branchInformationSize; b++) {
            segmentLength = Table.get("Branch length", b, "Branch information");
            segmentEuclidean = Table.get("Euclidean distance", b, "Branch information");
            segmentIntensity = Table.get("average intensity", b, "Branch information");

            if (segmentLength > minBranchLength &&
                segmentEuclidean != 0 &&
                segmentIntensity != 0) {
                nSkeletonSegments++;
            }
        }

        AnalyzeSkeletonResultsTableName = "Results";

        if (saveGroupsFilaments) {
            Table.save(groupedFilamentsOutputDir + "Grouped_filament_" + baseName + ".csv", AnalyzeSkeletonResultsTableName);
        }

        selectImage(skeletonAnalysisImageID);

		if (useMainPathDefinition)
    		helperArgument = "definition=main";
		else
    		helperArgument = "definition=segments";

		if (excludeSideBranchesFromOrientation)
    		helperArgument = helperArgument + ";exclude_branches=true";
		else
    		helperArgument = helperArgument + ";exclude_branches=false";

		if (saveIndividualFilaments) {
    		selectImage(skeletonAnalysisImageID);
    		run("CytoOrg Filament Bridge", helperArgument);

    		Table.save(individualFilamentsOutputDir + "Single_filament_" + baseName + ".csv", "CytoOrg Filaments v1.0");
    		close("CytoOrg Filaments v1.0");
    		close("CytoOrg Orientation v1.0");
		}

		selectImage(skeletonAnalysisImageID);
		run("CytoOrg Filament Bridge", helperArgument + ";filter_intensity=true");

		FilamentTableName = "CytoOrg Filaments v1.0";
		OrientationTableName = "CytoOrg Orientation v1.0";

		size = Table.size(FilamentTableName);

        nBranches = 0;
        sumLength = 0;
        sumLength2 = 0;
        maxLength = 0;
        sumEuclidean = 0;
        nShort = 0;
        nIntermediate = 0;
        nLong = 0;

        for (b = 0; b < size; b++) {
            Length = Table.get("Filament length", b, FilamentTableName);
            Euclidean = Table.get("Euclidean distance", b, FilamentTableName);

            if (Length > minBranchLength && Euclidean > 0) {
                nBranches++;
                sumLength += Length;
                sumLength2 += Length * Length;
                sumEuclidean += Euclidean;

                if (Length > maxLength)
                    maxLength = Length;

                if (StratifyLengths) {
                    if (Length < ShortIntermediate)
                        nShort++;
                    else if (Length < IntermediateLong)
                        nIntermediate++;
                    else
                        nLong++;
                }
            }
        }

        if (nBranches > 0) {
            meanLength = sumLength / nBranches;
            meanTortuosity = sumLength / sumEuclidean;
        }
        else {
            meanLength = NaN;
            meanTortuosity = NaN;
        }

        if (nBranches > 1) {
            variance = (sumLength2 - sumLength * sumLength / nBranches) / (nBranches - 1);

            if (variance < 0 && variance > -0.000000000001)
                variance = 0;

            stdLength = sqrt(variance);
        }
        else if (nBranches == 1) {
            stdLength = 0;
        }
        else {
            stdLength = NaN;
        }

        if (OrientationAnalysis) {
            sumOrientationCos2 = 0;
            sumOrientationSin2 = 0;
            sumOrientationCos4 = 0;
            sumOrientationSin4 = 0;
            orientationTotalWeight = 0;
            angularBinWeights = newArray(180 / orientationBinWidth);
            orientationTableSize = Table.size(OrientationTableName);

            for (b = 0; b < orientationTableSize; b++) {
                Length = Table.get("Filament length", b, OrientationTableName);
                Euclidean = Table.get("Euclidean distance", b, OrientationTableName);
                V1x = Table.get("V1 x", b, OrientationTableName);
                V1y = Table.get("V1 y", b, OrientationTableName);
                V2x = Table.get("V2 x", b, OrientationTableName);
                V2y = Table.get("V2 y", b, OrientationTableName);

                if (Length > minBranchLength && Euclidean > 0) {
                    filamentTortuosity = Length / Euclidean;

                    if (filamentTortuosity <= maxOrientationTortuosity) {
                        dx = V2x - V1x;
                        dy = V1y - V2y;
                        AngleWithXaxis = atan2(dy, dx);

                        if (weightOrientationByLength)
                            orientationWeight = Length;
                        else
                            orientationWeight = 1;

                        sumOrientationCos2 += orientationWeight * cos(2 * AngleWithXaxis);
                        sumOrientationSin2 += orientationWeight * sin(2 * AngleWithXaxis);
                        sumOrientationCos4 += orientationWeight * cos(4 * AngleWithXaxis);
                        sumOrientationSin4 += orientationWeight * sin(4 * AngleWithXaxis);
                        orientationTotalWeight += orientationWeight;

                        angleDegrees = AngleWithXaxis * 180 / PI;

                        while (angleDegrees < 0)
                            angleDegrees += 180;

                        while (angleDegrees >= 180)
                            angleDegrees -= 180;

                        angleBin = floor(angleDegrees / orientationBinWidth);
                        angularBinWeights[angleBin] += orientationWeight;
                    }
                }
            }

            if (orientationTotalWeight > 0) {
                meanCos2 = sumOrientationCos2 / orientationTotalWeight;
                meanSin2 = sumOrientationSin2 / orientationTotalWeight;
                meanCos4 = sumOrientationCos4 / orientationTotalWeight;
                meanSin4 = sumOrientationSin4 / orientationTotalWeight;

                orientationOrder = sqrt(meanCos2 * meanCos2 + meanSin2 * meanSin2);
                meanDoubledAngle = atan2(meanSin2, meanCos2);
                meanOrientation = 0.5 * meanDoubledAngle;
                meanOrientationDegrees = meanOrientation * 180 / PI;

                if (meanOrientationDegrees < 0)
                    meanOrientationDegrees += 180;

                circularSkewness = meanSin4 * cos(2 * meanDoubledAngle) -
                                   meanCos4 * sin(2 * meanDoubledAngle);
                circularKurtosis = meanCos4 * cos(2 * meanDoubledAngle) +
                                   meanSin4 * sin(2 * meanDoubledAngle);
            }
            else {
                orientationOrder = NaN;
                meanOrientationDegrees = NaN;
                circularSkewness = NaN;
                circularKurtosis = NaN;
            }

            if (saveAngularDistribution) {
                angularDistributionPath = angularDistributionOutputDir + "Angular_Distribution_" + baseName + ".tif";
                saveAngularDistributionPlot(angularBinWeights, orientationBinWidth,
                                            angularPlotMaximum, angularDistributionPath);
            }
        }

        close(FilamentTableName);
        close(OrientationTableName);
        close("Branch information");
        close(AnalyzeSkeletonResultsTableName);
    	
 	if(AnalyzeBranching){
 		
 		selectImage("Longest shortest paths");
    	
    	w = getWidth();
		h = getHeight();

		for (y = 1; y < h-1; y++) {
    		for (x = 1; x < w-1; x++) {
    
			if (getPixel(x,y) == 255) {
    			connected = false;

    		for (dy=-1; dy<=1; dy++) {
        		for (dx=-1; dx<=1; dx++) {

            if (dx==0 && dy==0)
                continue;

            if (getPixel(x+dx, y+dy) == 255)
                connected = true;
        	}
    	}

    		if (!connected)
        		setPixel(x, y, 96);
			}
   		}
	}
    	    	
        if (saveBranchMap || saveEnhancedBranchMap) {
			
			run("Duplicate...", "title=temp");
			run("RGB Color");
			
			for (y=0; y<h; y++) {
  			  for (x=0; x<w; x++) {

        		value = getPixel(x,y);

        		r = (value & 0xff0000) >> 16;
      			g = (value & 0x00ff00) >> 8;
       			b = value & 0x0000ff;

	        if (r==255 && g==255 && b==255) {
    	        setPixel(x,y,0xffff00);
        		}

	        if (r==195 && g==0 && b==93) {
    	        setPixel(x,y,0x00ffff);
        		}
    		}
		}
		
		if (saveBranchMap) {
			saveAs("Tiff", branchingOutputDir + "Branching_" + baseName + ".tif");
			}	
		
		if (saveEnhancedBranchMap) {
			
			run("Maximum...", "radius=3");	
			
			for (y=0; y<h; y++) {
  			  for (x=0; x<w; x++) {

        		value2 = getPixel(x,y);

        		r = (value2 & 0xff0000) >> 16;
      			g = (value2 & 0x00ff00) >> 8;
       			b = value2 & 0x0000ff;

	        if (r==255 && g==255 && b==255) {
    	        setPixel(x,y,0x00ffff);
        		}
    		}
		}
				
			saveAs("Tiff", improvedBranchingOutputDir + "Improved_Branching_" + baseName + ".tif");
			}	
		
		close();	
		
		}	
		
		MainBranches = 0;
		SideBranches = 0;

		for (y = 1; y < h-1; y++) {
    		for (x = 1; x < w-1; x++) {

				if (getPixel(x,y) == 96) {
					MainBranches++;
				}

				if (getPixel(x,y) == 255) {
					SideBranches++;
				}
    		}
		}	

		run("Duplicate...", "title=temp");
		setThreshold(1, 255);
 		run("Analyze Particles...", "summarize");
 		MainBranchesCount = Table.get("Count", 0, "Summary");
		close();
	
		setThreshold(96, 96);
		run("Create Mask");
		run("Analyze Particles...", "summarize");
		MainBranchesDiveded = Table.get("Count", 1, "Summary");	
		close();
	
		SideBranchesDiveded = nSkeletonSegments - MainBranchesDiveded;
	
		setThreshold(255, 255);
		run("Create Mask");
		run("Analyze Particles...", "summarize");
		SideBranchesCount = Table.get("Count", 2, "Summary");

    	close("*");
    	close("Summary");
    	
 	}
 	
 	close("*");
 
    if(OrientationAnalysis){
    	
    	open(SkeletonPath);
    	
    	run("Analyze Particles...", "size=3-Infinity pixel clear summarize");
		count_total = Table.get("Count", 0, "Summary");

		gradient = -1;

		for (g = 0; g < methods.length; g++) {
   	 		if (gradientMethod == methods[g]) {
        		gradient = g;
        		break;
    		}
		}

		run("OrientationJ Analysis", "tensor=" + orientationSigma + " gradient=" + gradient + " color-survey=on hsb=on hue=Orientation sat=Constant bri=Original-Image radian=on ");

		if (saveOrientationMap) {
			saveAs("Tiff", orientationOutputDir + "Orientation_" + baseName + ".tif");
		}	
		
		if (saveEnhancedOrientationMap) {
			run("Duplicate...", "title=Temp");
			run("Maximum...", "radius=3");	
			saveAs("Tiff", improvedOrientationOutputDir + "Improved_Orientation_" + baseName + ".tif");
			close("Temp");
		}	

		run("HSB Stack");
		run("Stack to Images");

		count_segment = 0;
		counts = newArray(256);

		for (b = 0; b < 256; b++) {

			AngleInterval = floor(directionChangeThreshold*256/180);

    		min = b;
    		max = (b + AngleInterval - 1) % 256;
    		
			selectImage("Hue");
			run("Duplicate...", "title=Hue_temp");
			selectImage("Hue_temp");

    	if (min <= max) {

			if (b == 0){
			        
			setThreshold(0.1, max);
        	run("Convert to Mask");
        
        	run("Clear Results");
			run("Analyze Particles...", "size=3-Infinity pixel clear summarize");
			count_segment = Table.get("Count", b + 1, "Summary");
		}

		else {

        	setThreshold(min, max);
        	run("Convert to Mask");
        
        	run("Clear Results");
			run("Analyze Particles...", "size=3-Infinity pixel clear summarize");
			count_segment = Table.get("Count", b + 1, "Summary");
		}

    } else {

        	setThreshold(min, 255);
        	run("Convert to Mask");
        	rename("part1");

			selectImage("Hue");
			run("Duplicate...", "title=part2");
			selectImage("part2");

        	setThreshold(0.1, max);
        	run("Convert to Mask");

        	imageCalculator("OR create", "part1", "part2");
        
        	run("Clear Results");
			run("Analyze Particles...", "size=3-Infinity pixel clear summarize");
			count_segment = Table.get("Count", b + 1, "Summary");

        	close("part1");
        	close("part2");

    	}

		counts[b] = count_segment;
	
    	close();
    	
	}

	close("Summary");

	sumOrientation = newArray(AngleInterval);
	nBins = floor(256 / AngleInterval);

	for (shift = 0; shift < AngleInterval; shift++) {

    	totalOrientation = 0;

    	for (q = 0; q < nBins; q++) {

        	index = shift + q * AngleInterval;

        	totalOrientation += counts[index];

    	}

    	sumOrientation[shift] = totalOrientation;
		
	}

	meanNofSegments = 0;

	for (m = 0; m < AngleInterval; m++) {
    	meanNofSegments += sumOrientation[m];
	}

	meanNofSegments = meanNofSegments / AngleInterval;
    NofSegments = meanNofSegments - count_total;
    
  }
	
  close("*");

    CellCytoskeletonArea = CytoskeletonArea / CellArea * 100;
    
    JunctionOnBranches = totalJunctions/nBranches;
    
    FilamentDensity = nBranches/CellArea * 100;
    
    if(StratifyLengths){
    	ShortFilamentDensity = nShort/CellArea * 100;
    	IntermediateFilamentDensity = nIntermediate/CellArea * 100;
    	LongFilamentDensity = nLong/CellArea * 100;
    }
    
    if(AnalyzeBranching){
    	
    	Branching = (SideBranchesCount / (SideBranchesCount + MainBranchesCount))*100;
    
    	BranchFilamentLengthRatio = SideBranches/(SideBranches+MainBranches);
    
    	RebranchingRate = SideBranchesDiveded / SideBranchesCount;
    
    }
    
    if(OrientationAnalysis){
        if (sumLength > 0)
            N_change_in_direction = NofSegments/sumLength;
        else
            N_change_in_direction = NaN;
    }

	AreaHeader = "Cell Area (" + micro + "m" + sup2 + ")";
    MeanThicknessHeader = "Mean Thickness (" + micro + "m" + ")";
    Q1ThicknessHeader = "Q1 Thickness (" + micro + "m" + ")";
    MedianThicknessHeader = "Median Thickness (" + micro + "m" + ")";
    Q3ThicknessHeader = "Q3 Thickness (" + micro + "m" + ")";
    CytoskeletonAreaHeader = "Cytoskeleton Area (" + micro + "m" + sup2 + ")";
    AcskAcellHeader = "Acsk/Acell (%)";
    LengthHeader = "Mean filament length (" + micro + "m)";
	StdLengthHeader = "Std Dev filament length (" + micro + "m)";
	TotalLengthHeader = "Total filament length (" + micro + "m)";
	MaxLengthHeader = "Maximum filament length (" + micro + "m)";
	FilamentDensityHeader = "Filament Density per 100 " + micro + "m" + sup2;
    ShortHeader = "N. of Short Filaments (<" + ShortIntermediate + " " + micro + "m)";
    IntermediateHeader = "N. of Intermediate Filaments (" + ShortIntermediate + " - " + IntermediateLong + " " + micro + "m)"; 
    LongHeader = "N. of Long Filaments (>" + IntermediateLong + " " + micro + "m)";
	ShortDensityHeader = "Short Filament per 100 " + micro + "m" + sup2;
	IntermediateDensityHeader = "Intermediate Filament per 100 " + micro + "m" + sup2;
	LongDensityHeader = "Long Filaments per 100 " + micro + "m" + sup2;
	DirectionChangeHeader = "Direction Changes per " + micro + "m"; 
    
    Table.set("Image", resultRow, baseName, "Cytoskeletal Organization - Results");

    Table.set(AreaHeader, resultRow, CellArea, "Cytoskeletal Organization - Results");
	Table.set("Cell Roundness", resultRow, Roundness, "Cytoskeletal Organization - Results");
	Table.set("Cell Solidity", resultRow, Solidity, "Cytoskeletal Organization - Results");	
    
	if(measureThickness){
		
		Table.set(MeanThicknessHeader, resultRow, meanThickness_um, "Cytoskeletal Organization - Results");
		Table.set(Q1ThicknessHeader, resultRow, q1_um, "Cytoskeletal Organization - Results");
		Table.set(MedianThicknessHeader, resultRow, median_um, "Cytoskeletal Organization - Results");
		Table.set(Q3ThicknessHeader, resultRow, q3_um, "Cytoskeletal Organization - Results");
	
	}

	Table.set(CytoskeletonAreaHeader, resultRow, CytoskeletonArea, "Cytoskeletal Organization - Results");
	Table.set(AcskAcellHeader, resultRow, CellCytoskeletonArea, "Cytoskeletal Organization - Results");
	
	Table.set("N. of Filaments", resultRow, nBranches, "Cytoskeletal Organization - Results");
	Table.set("N. of Junctions", resultRow, totalJunctions, "Cytoskeletal Organization - Results");
	Table.set("Junction/Filaments", resultRow, JunctionOnBranches, "Cytoskeletal Organization - Results");
	Table.set(LengthHeader, resultRow, meanLength, "Cytoskeletal Organization - Results");
	Table.set(StdLengthHeader, resultRow, stdLength, "Cytoskeletal Organization - Results");
	Table.set(TotalLengthHeader, resultRow, sumLength, "Cytoskeletal Organization - Results");
	Table.set(MaxLengthHeader, resultRow, maxLength, "Cytoskeletal Organization - Results");
	Table.set("Mean Tortuosity", resultRow, meanTortuosity, "Cytoskeletal Organization - Results");

	Table.set(FilamentDensityHeader, resultRow, FilamentDensity, "Cytoskeletal Organization - Results");
	
	if(StratifyLengths){
	Table.set(ShortHeader, resultRow, nShort, "Cytoskeletal Organization - Results");
	Table.set(IntermediateHeader, resultRow, nIntermediate, "Cytoskeletal Organization - Results");
	Table.set(LongHeader, resultRow, nLong, "Cytoskeletal Organization - Results");
	Table.set(ShortDensityHeader, resultRow, ShortFilamentDensity, "Cytoskeletal Organization - Results");
	Table.set(IntermediateDensityHeader, resultRow, IntermediateFilamentDensity, "Cytoskeletal Organization - Results");
	Table.set(LongDensityHeader, resultRow, LongFilamentDensity, "Cytoskeletal Organization - Results");
	}
	
	if(AnalyzeBranching){
		
		Table.set("Fraction of Branched Filaments (%)", resultRow, Branching, "Cytoskeletal Organization - Results");
		Table.set("Branches/Filaments length ratio", resultRow, BranchFilamentLengthRatio, "Cytoskeletal Organization - Results");
		Table.set("Rebranching rate", resultRow, RebranchingRate, "Cytoskeletal Organization - Results");
	
	}

	if(OrientationAnalysis){
		Table.set("Mean Orientation (Degrees)", resultRow, meanOrientationDegrees, "Cytoskeletal Organization - Results");
		Table.set("Orientation Order Parameter", resultRow, orientationOrder, "Cytoskeletal Organization - Results");
		Table.set("Circular Skewness", resultRow, circularSkewness, "Cytoskeletal Organization - Results");
		Table.set("Circular Kurtosis", resultRow, circularKurtosis, "Cytoskeletal Organization - Results");
		Table.set(DirectionChangeHeader, resultRow, N_change_in_direction, "Cytoskeletal Organization - Results");	
	}

	Table.update("Cytoskeletal Organization - Results");
	
	resultRow++;

}

if (MultiCellAnalysis && temporaryDir != "") {
    deleteTemporaryFolder(temporaryDir);
}

resultsPath = outputDir + "Results.csv";

Table.save(resultsPath, "Cytoskeletal Organization - Results");

csvContent = File.openAsString(resultsPath);
File.saveString(fromCharCode(65279) + csvContent, resultsPath);

setBatchMode(false);
restoreSettings;
print("Analysis Completed");

function getDimensionMismatchMessage(originalPath, maskPath, imageName) {

    open(originalPath);
    getDimensions(originalWidth, originalHeight, originalChannels, originalSlices, originalFrames);
    close();

    open(maskPath);
    getDimensions(maskWidth, maskHeight, maskChannels, maskSlices, maskFrames);
    close();

    if (originalWidth != maskWidth || originalHeight != maskHeight) {
        return "Image and mask dimensions do not match:\\n" +
               imageName + " = " + originalWidth + " x " + originalHeight + "\\n" +
               "segcell_" + imageName + " = " + maskWidth + " x " + maskHeight;
    }

    return "";
}

function createRoiControlImage(originalPath, baseName, temporaryDir, nCells) {

    controlFileName = "ROI_overlay_" + baseName + ".tif";
    controlFilePath = temporaryDir + controlFileName;

    open(originalPath);
    originalImageID = getImageID();

    Overlay.remove;

    for (controlIndex = 0; controlIndex < nCells; controlIndex++) {
        roiManager("Select", controlIndex);
        Overlay.addSelection("yellow", 2);
    }

    run("Select None");
    Overlay.useNamesAsLabels(false);
    Overlay.drawLabels(true);
    Overlay.setLabelColor("yellow");
    Overlay.setLabelFontSize(14, "bold back");
    Overlay.flatten;

    saveAs("Tiff", controlFilePath);
    close();

    if (isOpen(originalImageID)) {
        selectImage(originalImageID);
        close();
    }

    return controlFileName;
}

function getCategoryOutputDirectory(outputDir, categoryName) {

    categoryDir = outputDir + categoryName + File.separator;

    if (!File.exists(categoryDir)) {
        File.makeDirectory(categoryDir);
    }

    return categoryDir;
}

function copyRoiControlImages(controlFiles, temporaryDir, destinationDir) {

    for (controlIndex = 0; controlIndex < controlFiles.length; controlIndex++) {

        sourcePath = temporaryDir + controlFiles[controlIndex];
        destinationPath = destinationDir + controlFiles[controlIndex];

        if (File.exists(destinationPath)) {
            fileDeleteStatus = File.delete(destinationPath);
        }

        File.copy(sourcePath, destinationPath);

        if (!File.exists(destinationPath)) {
            return "Could not save ROI control image:\\n" + destinationPath;
        }
    }

    return "";
}

function isSupportedImageFile(fileName) {

    lowerFileName = toLowerCase(fileName);

    if (endsWith(lowerFileName, ".tif")) return true;
    if (endsWith(lowerFileName, ".tiff")) return true;
    if (endsWith(lowerFileName, ".png")) return true;
    if (endsWith(lowerFileName, ".jpg")) return true;
    if (endsWith(lowerFileName, ".jpeg")) return true;

    return false;
}

function getBaseName(fileName) {

    base = fileName;
    extensionDot = lastIndexOf(base, ".");

    if (extensionDot != -1) {
        base = substring(base, 0, extensionDot);
    }

    return base;
}

function deleteTemporaryFolder(folder) {

    if (folder != "" && File.exists(folder)) {

        temporaryFiles = getFileList(folder);

        for (deleteIndex = 0; deleteIndex < temporaryFiles.length; deleteIndex++) {
            fileDeleteStatus = File.delete(folder + temporaryFiles[deleteIndex]);
        }

        fileDeleteStatus = File.delete(folder);
    }
}

function saveAngularDistributionPlot(weights, binWidth, fixedMaximumFrequency, outputPath) {

    plotSize = 900;
    plotCenter = plotSize / 2;
    maximumRadius = 350;
    angleLabelRadius = 395;
    nBins = weights.length;
    maximumWeight = 0;
    totalWeight = 0;
    scaleMaximumFrequency = 0;

    for (binIndex = 0; binIndex < nBins; binIndex++) {
        totalWeight += weights[binIndex];

        if (weights[binIndex] > maximumWeight) {
            maximumWeight = weights[binIndex];
        }
    }

    if (totalWeight > 0) {
        if (fixedMaximumFrequency > 0)
            scaleMaximumFrequency = fixedMaximumFrequency;
        else
            scaleMaximumFrequency = 100 * maximumWeight / totalWeight;
    }

    newImage("Angular Distribution", "RGB white", plotSize, plotSize, 1);
    autoUpdate(false);
    setLineWidth(1);

    for (ring = 1; ring <= 4; ring++) {
        ringRadius = maximumRadius * ring / 4;
        setColor("#D7D7D7");
        drawOval(round(plotCenter - ringRadius), round(plotCenter - ringRadius),
                 round(2 * ringRadius), round(2 * ringRadius));
    }

    setColor("#E3E3E3");

    for (angleDegrees = 0; angleDegrees < 360; angleDegrees += 30) {
        angleRadians = angleDegrees * PI / 180;
        spokeX = plotCenter + maximumRadius * cos(angleRadians);
        spokeY = plotCenter - maximumRadius * sin(angleRadians);
        drawLine(plotCenter, plotCenter, round(spokeX), round(spokeY));
    }

    setColor("#7A7A7A");
    setLineWidth(2);
    drawOval(plotCenter - maximumRadius, plotCenter - maximumRadius,
             2 * maximumRadius, 2 * maximumRadius);

    if (maximumWeight > 0 && scaleMaximumFrequency > 0) {
        setColor("#1F4E79");
        setLineWidth(4);
        previousX = 0;
        previousY = 0;

        for (pointIndex = 0; pointIndex <= 2 * nBins; pointIndex++) {
            binIndex = pointIndex % nBins;
            turnOffset = floor(pointIndex / nBins) * 180;

            if (pointIndex == 2 * nBins) {
                turnOffset = 360;
            }

            angleDegrees = (binIndex + 0.5) * binWidth + turnOffset;
            angleRadians = angleDegrees * PI / 180;
            binFrequency = 100 * weights[binIndex] / totalWeight;
            radialDistance = maximumRadius * binFrequency / scaleMaximumFrequency;

            if (radialDistance > maximumRadius)
                radialDistance = maximumRadius;

            pointX = plotCenter + radialDistance * cos(angleRadians);
            pointY = plotCenter - radialDistance * sin(angleRadians);

            if (pointIndex > 0) {
                drawLine(round(previousX), round(previousY), round(pointX), round(pointY));
            }

            previousX = pointX;
            previousY = pointY;
        }
    }

    setColor("#333333");
    setFont("SansSerif", 18, "antialiased");
    setJustification("center");
    degreeSymbol = fromCharCode(176);

    for (angleDegrees = 0; angleDegrees < 360; angleDegrees += 30) {
        angleRadians = angleDegrees * PI / 180;
        labelX = plotCenter + angleLabelRadius * cos(angleRadians);
        labelY = plotCenter - angleLabelRadius * sin(angleRadians) + 7;
        drawString(angleDegrees + degreeSymbol, round(labelX), round(labelY));
    }

    if (totalWeight > 0) {
        setFont("SansSerif", 16, "antialiased");
        setJustification("left");

        for (ring = 1; ring <= 4; ring++) {
            ringRadius = maximumRadius * ring / 4;
            ringFrequency = scaleMaximumFrequency * ring / 4;
            drawString(d2s(ringFrequency, 1) + "%", plotCenter + 8, round(plotCenter - ringRadius + 6));
        }
    }

    updateDisplay();
    autoUpdate(true);
    saveAs("Tiff", outputPath);
    close();
}

'''
    IJ.runMacro(macroSource)
}
finally {
    Executer.removeCommandListener(bridge)
}