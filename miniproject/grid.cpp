#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("membrane");

  double a = 1;
  
  double lc = 1e-2;
  gmsh::model::geo::addPoint(0,0,0, lc, 0);
  gmsh::model::geo::addPoint(a,0,0, lc, 1);
  gmsh::model::geo::addPoint(a,a,0, lc, 2);
  gmsh::model::geo::addPoint(0,a,0, lc, 3);

  gmsh::model::geo::addLine(0, 1, 0);
  gmsh::model::geo::addLine(1, 2, 1);
  gmsh::model::geo::addLine(2, 3, 2);
  gmsh::model::geo::addLine(3, 0, 3);


  gmsh::model::geo::addCurveLoop({0, 1, 2, 3}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);


  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(3);

  gmsh::write("membrane.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}

