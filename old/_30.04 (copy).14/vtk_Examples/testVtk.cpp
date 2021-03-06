#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

int main()
{
    int nx = 10, ny = 10, nz = 10;

    vtkSmartPointer<vtkImageData> imageData =
            vtkSmartPointer<vtkImageData>::New();

    imageData->SetDimensions(nx, ny, nz);

    vtkSmartPointer<vtkDoubleArray> director =
            vtkSmartPointer<vtkDoubleArray>::New();

    director->SetNumberOfComponents(3);
    director->SetNumberOfTuples(nx * ny * nz);

    vtkSmartPointer<vtkDoubleArray> energy =
            vtkSmartPointer<vtkDoubleArray>::New();

    energy->SetNumberOfComponents(1);
    energy->SetNumberOfTuples(nx * ny * nz);

    for (int i = 0; i < director->GetNumberOfTuples(); ++i) {
        double t = 1.0;
        double p = 0.0;
        double e = 5.0;
        double x = sin(t) * cos(p),
                y = sin(t) * sin(p),
                z = cos(t);

        director->SetTuple3(i, x, y, z);
        energy->SetValue(i, e);
    }

    imageData->GetPointData()->AddArray(director);
    director->SetName("Director");

    imageData->GetPointData()->AddArray(energy);
    energy->SetName("Energy");

    vtkSmartPointer<vtkXMLImageDataWriter> writer =
            vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetFileName("test.vti");
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(imageData->GetProducerPort());
#else
    writer->SetInputData(imageData);
#endif
    writer->Write();

    return EXIT_SUCCESS;
}
