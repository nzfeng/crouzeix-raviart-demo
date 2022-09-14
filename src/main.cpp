#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
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

float TCOEF = 1.0;

/*
 * Convert edge midpoint scalar data to vertex data by averaging.
 */
Vector<double> edgeMidpointDataToVertexData(const Vector<double>& u) {

  geometry->requireVertexIndices();
  geometry->requireEdgeIndices();
  Vector<double> w = Vector<double>::Zero(mesh->nVertices());
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
  // Pick two random point sources.
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> distr(0, E - 1); // distribution in range [0, E-1]
  size_t ePosIdx = distr(rng);
  size_t eNegIdx = distr(rng);
  Edge ePos = mesh->edge(ePosIdx);
  Edge eNeg = mesh->edge(eNegIdx);
  Vector<double> rho = Vector<double>::Zero(E);
  rho[ePos] = 1;
  rho[eNeg] = -1;
  std::cerr << "ePos: " << ePos << std::endl;
  std::cerr << "eNeg: " << eNeg << std::endl;
  // Need the r.h.s. to integrate to zero.
  geometry->requireFaceAreas();
  double totalArea = 0.;
  for (Face f : mesh->faces()) {
    totalArea += geometry->faceAreas[f];
  }
  geometry->unrequireFaceAreas();
  Vector<double> rhoBar = Vector<double>::Ones(E) * (M * rho).sum() / totalArea;
  Vector<double> RHS = M * (rhoBar - rho);

  Vector<double> sol = solvePositiveDefinite(L, RHS);
  // Display on vertices
  psMesh->addVertexScalarQuantity("u", edgeMidpointDataToVertexData(sol));

  geometry->unrequireCrouzeixRaviartLaplacian();
  geometry->unrequireCrouzeixRaviartMassMatrix();
}

/*
 * Diffuse vectors using the Crouzeix-Raviart connection Laplacian.
 */
void diffuseVectors() {

  geometry->requireCrouzeixRaviartConnectionLaplacian();

  SparseMatrix<std::complex<double>> L = geometry->crouzeixRaviartConnectionLaplacian;

  // TODO: Pick a curve and place vectors along it

  geometry->unrequireCrouzeixRaviartConnectionLaplacian();
}

void myCallback() {

  ImGui::SliderFloat("tCoef", &TCOEF, 0., 100.);

  if (ImGui::Button("Solve Poisson problem")) {
    solvePoissonProblem();
  }

  if (ImGui::Button("Diffuse vectors")) {
    diffuseVectors();
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

  // TODO: Set up solver

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}