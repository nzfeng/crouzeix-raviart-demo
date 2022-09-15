#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <random>

#include "args/args.hxx"
#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// A bunch of random points are selected upon program startup.
int nSOURCE_POINTS = 2;
std::vector<size_t> SOURCE_POINTS;
std::vector<double> SOURCE_STRENGTHS;
// A bunch of random edges are selected upon program startup; unit vectors placed at the midpoints of these edges are
// used as the source for the diffuseVectors() function.
std::vector<Edge> SELECTED_EDGES;
std::vector<double> SOURCE_MAGNITUDES;
std::vector<double> SOURCE_ANGLES;
int nEDGES_SELECTED = 1;
double MEAN_EDGE_LENGTH = 0.;


/*
 * Convert edge midpoint scalar data to vertex data by averaging.
 */
template <typename T>
Vector<T> edgeMidpointDataToVertexData(const Vector<T>& u) {

  geometry->requireVertexIndices();
  geometry->requireEdgeIndices();
  Vector<T> w = Vector<T>::Zero(mesh->nVertices());
  for (Vertex v : mesh->vertices()) {
    size_t vIdx = geometry->vertexIndices[v];
    for (Edge e : v.adjacentEdges()) {
      size_t eIdx = geometry->edgeIndices[e];
      w[vIdx] += u[eIdx];
    }
    w[vIdx] /= v.degree();
  }
  geometry->requireVertexIndices();
  geometry->unrequireEdgeIndices();

  return w;
}

// /*
//  * Need to average vectors onto vertices by first converting to the vertex tangent space.
//  */
// Vector<std::complex<double>> edgeMidpointDataToVertexData(const Vector<std::complex<double>>& u) {

//   geometry->requireVertexIndices();
//   geometry->requireEdgeIndices();
//   Vector<std::complex<double>> w = Vector<std::complex<double>>::Zero(mesh->nVertices());
//   for (Vertex v : mesh->vertices()) {
//     size_t vIdx = geometry->vertexIndices[v];
//     for (Edge e : v.adjacentEdges()) {
//       size_t eIdx = geometry->edgeIndices[e];
//       w[vIdx] += u[eIdx];
//     }
//     w[vIdx] /= v.degree();
//   }
//   geometry->requireVertexIndices();
//   geometry->unrequireEdgeIndices();

//   return w;
// }

void displaySourcePoints() {

  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  for (size_t i = 0; i < nSOURCE_POINTS; i++) {
    Edge e = mesh->edge(SOURCE_POINTS[i]);
    Vector3 a = geometry->vertexPositions[e.firstVertex()];
    Vector3 b = geometry->vertexPositions[e.secondVertex()];
    positions.push_back(0.5 * (a + b));
  }
  psMesh->addSurfaceGraphQuantity("source points", positions, edgeInds)->setColor({0.0, 0.0, 0.0})->setEnabled(true);
}

/*
 * Solve a Poisson equation with Dirichlet sources using the Crouzeix-Raviart Laplacian.
 */
void solvePoissonProblem() {

  geometry->requireCrouzeixRaviartLaplacian();
  geometry->requireCrouzeixRaviartMassMatrix();

  SparseMatrix<double> L = geometry->crouzeixRaviartLaplacian;
  SparseMatrix<double> M = geometry->crouzeixRaviartMassMatrix;

  size_t E = mesh->nEdges();
  // Set source vector.
  Vector<double> rho = Vector<double>::Zero(E);
  for (size_t i = 0; i < nSOURCE_POINTS; i++) {
    rho[SOURCE_POINTS[i]] = SOURCE_STRENGTHS[i];
    std::cerr << "Edge " << SOURCE_POINTS[i] << "\tStrength: " << SOURCE_STRENGTHS[i] << std::endl;
  }
  // Need the r.h.s. to integrate to zero.
  geometry->requireFaceAreas();
  double totalArea = 0.;
  for (Face f : mesh->faces()) {
    totalArea += geometry->faceAreas[f];
  }
  geometry->unrequireFaceAreas();
  Vector<double> rhoBar = Vector<double>::Ones(E) * (M * rho).sum() / totalArea;
  Vector<double> RHS = M * (rho - rhoBar);

  Vector<double> sol = solvePositiveDefinite(L, RHS);
  // Display on vertices
  size_t V = mesh->nVertices();
  Vector<double> sol_verts = edgeMidpointDataToVertexData(sol);
  sol_verts -= Vector<double>::Ones(mesh->nVertices()) * sol_verts.minCoeff();
  psMesh->addVertexScalarQuantity("u", sol_verts)->setEnabled(true);
  displaySourcePoints();

  geometry->unrequireCrouzeixRaviartLaplacian();
  geometry->unrequireCrouzeixRaviartMassMatrix();
}

void displaySourceVectors() {

  // Display the origin points
  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  for (Edge e : SELECTED_EDGES) {
    Vector3 a = geometry->vertexPositions[e.firstVertex()];
    Vector3 b = geometry->vertexPositions[e.secondVertex()];
    Vector3 midpoint = 0.5 * (a + b);
    positions.push_back(midpoint);
  }
  auto points = psMesh->addSurfaceGraphQuantity("vector origins", positions, edgeInds);
  points->setColor({0.0, 0.0, 0.0});
  points->setEnabled(true);

  // Display the direction by drawing a segment in the direction
  TraceGeodesicResult tracedGeodesic;
  TraceOptions traceOptions;
  traceOptions.includePath = true;
  std::vector<std::array<size_t, 2>> edgeInds_vectors;
  std::vector<Vector3> positions_vectors;
  for (size_t i = 0; i < nEDGES_SELECTED; i++) {
    Vector2 vec = SOURCE_MAGNITUDES[i] * MEAN_EDGE_LENGTH / 2 * Vector2::fromAngle(SOURCE_ANGLES[i]);
    tracedGeodesic = traceGeodesic(*geometry, SurfacePoint(SELECTED_EDGES[i], 0.5), vec, traceOptions);
    std::vector<SurfacePoint>& pathPoints = tracedGeodesic.pathPoints;
    int offset = positions_vectors.size();
    for (SurfacePoint& pt : pathPoints) {
      positions_vectors.push_back(pt.interpolate(geometry->vertexPositions));
    }
    for (size_t j = 0; j < pathPoints.size() - 1; j++) {
      edgeInds_vectors.push_back({offset + j, offset + j + 1});
    }
  }
  auto vectors = psMesh->addSurfaceGraphQuantity("vectors", positions_vectors, edgeInds_vectors);
  vectors->setColor({1.0, 0.0, 1.0});
  vectors->setEnabled(true);
}

void displayTransportedVectors(const Vector<std::complex<double>>& sol) {

  TraceGeodesicResult tracedGeodesic;
  TraceOptions traceOptions;
  traceOptions.includePath = true;
  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  size_t E = mesh->nEdges();
  for (size_t i = 0; i < E; i++) {
    Vector2 vec = Vector2::fromComplex(sol[i]);
    tracedGeodesic =
        traceGeodesic(*geometry, SurfacePoint(mesh->edge(i), 0.5), MEAN_EDGE_LENGTH / 2. * vec, traceOptions);
    std::vector<SurfacePoint>& pathPoints = tracedGeodesic.pathPoints;
    int offset = positions.size();
    for (SurfacePoint& pt : pathPoints) {
      positions.push_back(pt.interpolate(geometry->vertexPositions));
    }
    for (size_t j = 0; j < pathPoints.size() - 1; j++) {
      edgeInds.push_back({offset + j, offset + j + 1});
    }
  }
  auto vectors = psMesh->addSurfaceGraphQuantity("transported vectors", positions, edgeInds);
  vectors->setColor({0.0, 1.0, 0.0});
  vectors->setEnabled(true);
}

/*
 * Do the Vector Heat Method.
 */
void doParallelTransport() {

  geometry->requireCrouzeixRaviartLaplacian();
  geometry->requireCrouzeixRaviartConnectionLaplacian();
  geometry->requireCrouzeixRaviartMassMatrix();
  geometry->requireEdgeIndices();

  size_t E = mesh->nEdges();
  size_t V = mesh->nVertices();

  SparseMatrix<std::complex<double>>& L = geometry->crouzeixRaviartConnectionLaplacian;
  SparseMatrix<double>& M = geometry->crouzeixRaviartMassMatrix;
  double shortTime = MEAN_EDGE_LENGTH * MEAN_EDGE_LENGTH;
  SparseMatrix<std::complex<double>> LHS = M.cast<std::complex<double>>() + shortTime * L; // L is positive

  SparseMatrix<double> LHS_scalar = M + shortTime * geometry->crouzeixRaviartLaplacian;
  PositiveDefiniteSolver<double> scalarHeatSolver(LHS_scalar);
  Vector<double> dataRHS = Vector<double>::Zero(E);
  Vector<double> indicatorRHS = Vector<double>::Zero(E);

  // Set up r.h.s.
  Vector<std::complex<double>> rhs_vec = Vector<std::complex<double>>::Zero(E);
  for (size_t i = 0; i < nEDGES_SELECTED; i++) {
    size_t eIdx = geometry->edgeIndices[SELECTED_EDGES[i]];
    std::complex<double> vec = Vector2::fromAngle(SOURCE_ANGLES[i]);
    rhs_vec[eIdx] = vec;
    dataRHS[eIdx] = SOURCE_MAGNITUDES[i];
    indicatorRHS[eIdx] = 1;
  }

  // The C-R connection Laplacian should always be PSD, since z^T L z takes the inner product of two PL vector fields
  // within a face, and inner products are PSD
  Vector<std::complex<double>> sol = solvePositiveDefinite(LHS, rhs_vec);

  // Get the right magnitudes.
  Vector<double> dataSol = scalarHeatSolver.solve(dataRHS);
  Vector<double> indicatorSol = scalarHeatSolver.solve(indicatorRHS);
  Vector<double> interpResult = dataSol.array() / indicatorSol.array();
  for (size_t i = 0; i < E; i++) {
    sol[i] = interpResult[i] * Vector2::fromComplex(sol[i]).normalize();
  }

  // Display initial vectors and solution.
  displaySourceVectors();

  // TODO: Currently no Polyscope function SurfaceMesh::addEdgeIntrinsicVectorQuantity() to visualize vectors at edges.
  // For now, do a hacky version.
  displayTransportedVectors(sol);

  geometry->unrequireCrouzeixRaviartLaplacian();
  geometry->unrequireCrouzeixRaviartConnectionLaplacian();
  geometry->unrequireCrouzeixRaviartMassMatrix();
  geometry->unrequireEdgeIndices();
}

void myCallback() {

  // ImGui::SliderFloat("tCoef", &TCOEF, 0., 100.);

  if (ImGui::Button("Solve Poisson problem")) {
    psMesh->removeAllQuantities();
    solvePoissonProblem();
  }

  if (ImGui::Button("Parallel transport")) {
    psMesh->removeAllQuantities();
    doParallelTransport();
  }
}

int main(int argc, char** argv) {

  // Configure the argument parser
  args::ArgumentParser parser("crouzeix-raviart-demo");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help& h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(args::get(inputFilename)),
                                          geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                          polyscopePermutations(*mesh));

  // Compute the mean edge length of the mesh.
  MEAN_EDGE_LENGTH = 0.;
  geometry->requireEdgeLengths();
  for (Edge e : mesh->edges()) {
    MEAN_EDGE_LENGTH += geometry->edgeLengths[e];
  }
  MEAN_EDGE_LENGTH /= mesh->nEdges();
  geometry->unrequireEdgeLengths();

  size_t E = mesh->nEdges();
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> distr(0, E - 1);
  std::uniform_real_distribution<float> distrStrength(-2, 2);
  // Pick a bunch of random points and their magnitudes.
  for (size_t i = 0; i < nSOURCE_POINTS; i++) {
    SOURCE_POINTS.push_back(distr(rng));
    SOURCE_STRENGTHS.push_back(distrStrength(rng));
  }
  // Pick a bunch of random edges to place source vectors on, pick a bunch of random magnitudes, and pick a bunch of
  // random angles.
  std::uniform_int_distribution<std::mt19937::result_type> distrMag(1, 4);
  std::uniform_real_distribution<float> distrAngle(0, 2. * PI);
  for (size_t i = 0; i < nEDGES_SELECTED; i++) {
    SELECTED_EDGES.push_back(mesh->edge(distr(rng)));
    SOURCE_MAGNITUDES.push_back(distrMag(rng));
    SOURCE_ANGLES.push_back(distrAngle(rng));
  }

  SELECTED_EDGES[0] = mesh->edge(9757);
  SOURCE_MAGNITUDES[0] = 2;
  SOURCE_ANGLES[0] = PI / 2.;

  // Set vertex tangent spaces
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
  }
  psMesh->setVertexTangentBasisX(vBasisX);

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}