#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("circle");

  double torus_radius = 2.0;
  double tube_radius = 0.5;
  
  double lc = 1e-1;
  gmsh::model::occ::addPoint(torus_radius, 0, 0, lc, 10);
  gmsh::model::occ::addPoint(torus_radius, tube_radius, 0, lc, 1);
  gmsh::model::occ::addPoint(torus_radius + tube_radius, 0, 0, lc, 2);
  gmsh::model::occ::addPoint(torus_radius, -tube_radius, 0, lc, 3);
  gmsh::model::occ::addPoint(torus_radius - tube_radius, 0, 0, lc, 4);
  
  gmsh::model::occ::addCircleArc(1, 10, 2, 1);
  gmsh::model::occ::addCircleArc(2, 10, 3, 2);
  gmsh::model::occ::addCircleArc(3, 10, 4, 3);
  gmsh::model::occ::addCircleArc(4, 10, 1, 4);

  gmsh::model::occ::addCurveLoop({1, 2, 3, 4}, 1);
  gmsh::model::occ::addPlaneSurface({1}, 1);

  gmsh::model::occ::synchronize();

  std::vector<std::pair<int, int>> out;
  gmsh::model::occ::revolve({{2, 1}},     // поверхность
                               0, 0, 0,         // точка на оси вращения
                               0, 1, 0,         // направление оси Y
                               2 * M_PI,        // полный оборот (360°)
                               out);
                               
                               
                               
  gmsh::model::occ::synchronize();

  gmsh::model::mesh::generate(3);

  gmsh::write("torus.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}

