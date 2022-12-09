#include <vtkPolyData.h>
#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>
#include <vtkCylinderSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
 
int main(int argc, char *argv[])
{
  if(argc != 2)
    {
    std::cout << "Required parameters: filename.stl" << std::endl;
    return EXIT_FAILURE;
    }
 
  std::string filename = argv[1];
 

  // Create a sphere
  vtkSmartPointer<vtkCylinderSource> cylinderSource =
    vtkSmartPointer<vtkCylinderSource>::New();
  cylinderSource->SetCenter(50.0, 50.0, 0.0);
  cylinderSource->SetRadius(6.5);
  cylinderSource->SetHeight(450.0);
  cylinderSource->SetResolution(1000);

 
 vtkSmartPointer<vtkCylinderSource> cylinderSource1 =
    vtkSmartPointer<vtkCylinderSource>::New();
  cylinderSource->SetCenter(150.0, 150.0, 0.0);
  cylinderSource->SetRadius(6.5);
  cylinderSource->SetHeight(450.0);
  cylinderSource->SetResolution(1000);

  vtkSmartPointer<vtkSTLWriter> stlWriter =
    vtkSmartPointer<vtkSTLWriter>::New();
  stlWriter->SetFileName(filename.c_str());
  stlWriter->SetInputConnection(cylinderSource->GetOutputPort());
stlWriter->SetInputConnection(cylinderSource1->GetOutputPort());
  stlWriter->Write();
 
  // Read and display for verification
  vtkSmartPointer<vtkSTLReader> reader =
    vtkSmartPointer<vtkSTLReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
 
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(cylinderSource->GetOutputPort());
  mapper->SetInputConnection(cylinderSource1->GetOutputPort());

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);


 
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  renderer->AddActor(actor);

  renderer->SetBackground(.3, .6, .3); // Background color green
 
  renderWindow->Render();
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
}
