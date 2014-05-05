/*=========================================================================
 
 Program: Visualization Toolkit
 Module: vtkXMLStructuredGridWriter.h
 
 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE. See the above copyright notice for more information.
 
=========================================================================*/
#ifndef __vtkXMLStructuredGridWriter_h
#define __vtkXMLStructuredGridWriter_h

#include "vtkIOXMLModule.h" // For export macro
#include "vtkXMLStructuredDataWriter.h"
 
class vtkStructuredGrid;
 
class VTKIOXML_EXPORT vtkXMLStructuredGridWriter : public vtkXMLStructuredDataWriter
{
public:
 static vtkXMLStructuredGridWriter* New();
 vtkTypeMacro(vtkXMLStructuredGridWriter,vtkXMLStructuredDataWriter);
 void PrintSelf(ostream& os, vtkIndent indent);

 //BTX

 vtkStructuredGrid* GetInput();
 //ETX

 const char* GetDefaultFileExtension();

protected:
 vtkXMLStructuredGridWriter();
 ~vtkXMLStructuredGridWriter();

 // see algorithm for more info
 virtual int FillInputPortInformation(int port, vtkInformation* info);

 void WriteAppendedPiece(int index, vtkIndent indent);
 void WriteAppendedPieceData(int index);
 void WriteInlinePiece(vtkIndent indent);
 void GetInputExtent(int* extent);
 const char* GetDataSetName();
 void CalculateSuperclassFraction(float* fractions);

 // The position of the appended data offset attribute for the points
 // array.
 OffsetsManagerGroup *PointsOM; //one per piece

 virtual void AllocatePositionArrays();
 virtual void DeletePositionArrays();

private:
 vtkXMLStructuredGridWriter(const vtkXMLStructuredGridWriter&); // Not implemented.
 void operator=(const vtkXMLStructuredGridWriter&); // Not implemented.
};

#endif
