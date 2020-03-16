import java.awt.*;
import Jama.*; 
import java.util.ArrayList;
import java.util.HashMap;


public class MyPolygon extends Polygon {
	// For triangulation 
	public Polygon[] triangles;
	public int[] tempX = new int[3];
	public int[] tempY = new int[3];
	public int numTriangles; 
	public boolean triangulated = false; 

	// Mesh 
	public int[][] initialVertices; 
	public int[][] deformedVertices; 
	public int numVertices;
	public ArrayList<Integer> constrainedIdx; 
	public double[][][] triangleLocal; // Vertices in triangle local coordinates

	// Matrix from step one
	public Matrix gPrimeInvB;
	HashMap<XYKey, Integer> vertMap; // Vertex Map (Compiled in step one)
	int[] vVector; // Vector for computing G 
	// Updated Polygon and Triangles 
	public Polygon[] updatedTri;
	public Polygon updatedPoly; 

	// public int[][] gComponentsSingle = new int[6][6];
	// public int[] constraintX = new int[3];
	// public int[] constraintY = new int[3]; 
	// public int numConstraints = 0; 
	// public int[][] gMatrix; 
	// public int[] v;
	// 
	// public int[] newPos; 

	private void resizeTriangles(int num) {
		assert(tempX.length == tempY.length); 
		while(3*num > tempX.length) {
			// Double array capacity 
			int[] copyX = new int[tempX.length*2]; 
			int[] copyY = new int[tempX.length*2];
			// Copy data across 
			System.arraycopy(tempX, 0, copyX, 0, tempX.length);
			System.arraycopy(tempY, 0, copyY, 0, tempX.length);
			// Reset references 
			tempX = copyX; 
			tempY = copyY; 
		}
		assert(tempX.length == tempY.length);
	}

	/*
	private void resizeConstraints(int num) {
		assert(constraintX.length == constraintY.length);
		while(num > constraintX.length) {
			// Double array capacity
			int[] copyOldX = new int[constraintX.length*2]; 
			int[] copyOldY = new int[constraintX.length*2];
			// Copy data across 
			System.arraycopy(constraintX, 0, copyOldX, 0, tempX.length);
			System.arraycopy(constraintY, 0, copyOldY, 0, tempX.length);
			
			// Reset references 
			constraintX = copyOldX;
			constraintY = copyOldY;
		}
		assert(constraintX.length == constraintY.length);
	}
	*/

	public MyPolygon() {}
	public MyPolygon(int x[], int y[], int n) {}

	// Left 
	public boolean left(int ax, int ay, int bx, int by, int cx, int cy) {
		return (bx - ax)*(cy - ay) - (cx - ax)*(by - ay) < 0;
	}

	// Dot 
	public int DotProd(int[] x, int[] y) {
		return (x[0]*y[0]) + (x[1]*y[1]);
	}

	// Square length 
	public double squareLength(int[] vec) {
		return ((double)vec[0]*vec[0]) + ((double)vec[1]*vec[1]);
	}

	// O(n^4) implementation 
	public void triangulate() {
		numTriangles = 0; 

		// Compute z = x^2 + y^2
		int[] z = new int[npoints]; 
		for (int i=0; i<npoints; i++) {
			z[i] = (xpoints[i]*xpoints[i]) +
				(ypoints[i]*ypoints[i]);
		}
		// For each triple (i,j,k)
		int xn, yn, zn;
		double cx, cy;
		boolean flag; 
		for (int i=0; i<npoints-2; i++) {
			for (int j=i+1; j<npoints; j++) {
				for (int k=i+1; k<npoints; k++) {
					// Compute normal 
					xn = (ypoints[j]-ypoints[i])*(z[k]-z[i]) -
						(ypoints[k]-ypoints[i])*(z[j]-z[i]);
					yn = (xpoints[k]-xpoints[i])*(z[j]-z[i]) -
						(xpoints[j]-xpoints[i])*(z[k]-z[i]);
					zn = (xpoints[j]-xpoints[i])*(ypoints[k]-ypoints[i]) -
						(xpoints[k]-xpoints[i])*(ypoints[j]-ypoints[i]);
					// Examine faces on the bottom of the paraboloid
					if(zn < 0) {
						flag = true; 
						for(int m=0; m<npoints; m++) {
							flag = flag &&
								((xpoints[m]-xpoints[i])*xn +
								(ypoints[m]-ypoints[i])*yn +
								(z[m]-z[i])*zn <= 0);
						}
						if (flag) {
							// Exclude triangles entirely outside 
							cx = (xpoints[i] + xpoints[j] + xpoints[k]) / 3.0;
							cy = (ypoints[i] + ypoints[j] + ypoints[k]) / 3.0;
							if(contains(cx,cy)) {
								numTriangles++; 
								resizeTriangles(numTriangles);
								tempX[3*(numTriangles-1)] = xpoints[i];
								tempX[3*(numTriangles-1) + 1] = xpoints[j];
								tempX[3*(numTriangles-1) + 2] = xpoints[k];
								tempY[3*(numTriangles-1)] = ypoints[i];
								tempY[3*(numTriangles-1) + 1] = ypoints[j];
								tempY[3*(numTriangles-1) + 2] = ypoints[k];
								tempY[3*(numTriangles-1) + 2] = ypoints[k];
							}
						}
					}
				}
			}
		}
		triangles = new Polygon[numTriangles];
		int[] newPolyX = new int[3]; 
		int[] newPolyY = new int[3];
		for(int i=0; i<numTriangles; i++) {
			for(int j=0; j<3; j++) {
				newPolyX[j] = tempX[3*i + j];
				newPolyY[j] = tempY[3*i + j]; 
			}
			triangles[i] = new Polygon(newPolyX, newPolyY, 3);
		}
		triangulated = true; 
	}

	// Initialize from Mesh
	public void initializeMesh() {
		constrainedIdx = new ArrayList<Integer>();
		constrainedIdx.clear();

		// Copy vertices 
		ArrayList<int[]> verts = new ArrayList<int[]>(); 
		boolean inVerts; 
		numVertices = 0; 
		for(int i=0; i<numTriangles; i++) { 
			Polygon triangle = triangles[i]; 
			for(int j=0; j<3; j++) {
				inVerts = false;
				for(int k=0; k<numVertices; k++) {
					if(triangle.xpoints[j] == verts.get(k)[0] && 
						triangle.ypoints[j] == verts.get(k)[1]) {
						inVerts = true; 
						break;
					}
				}
				if(!inVerts) {
					verts.add(new int[]{triangle.xpoints[j], triangle.ypoints[j]});
					numVertices++;
				}
			}
		}
		initialVertices = new int[numVertices][2];
		deformedVertices = new int[numVertices][2];
		for(int i=0; i<numVertices; i++) { 
			int[] aVertex = verts.get(i);
			initialVertices[i][0] = aVertex[0];
			initialVertices[i][1] = aVertex[1];
			deformedVertices[i][0] = aVertex[0];
			deformedVertices[i][1] = aVertex[1];
		}

		// Triangle in local coodinates 
		triangleLocal = new double[numTriangles][3][2]; 
		int n0, n1, n2; 
		int[] v0 = new int[2], v1 = new int[2], v2 = new int[2];
		int[] v01 = new int[2], v01Rot90 = new int[2];
		int[] v02 = new int[2];
		for(int i=0; i<numTriangles; i++) {
			for(int j=0; j<3; j++) {
				n0 = j; 
				n1 = (j+1)%3;
				n2 = (j+2)%3; 

				v0[0] = triangles[i].xpoints[n0];
				v0[1] = triangles[i].ypoints[n0];
				v1[0] = triangles[i].xpoints[n1];
				v1[1] = triangles[i].ypoints[n1];
				v2[0] = triangles[i].xpoints[n2];
				v2[1] = triangles[i].ypoints[n2];

				v01[0] = v1[0] - v0[0];
				v01[1] = v1[1] - v0[1];
				v01Rot90[0] = v01[1];
				v01Rot90[1] = -v01[0];

				v02[0] = v2[0] - v0[0];
				v02[1] = v2[1] - v0[1];

				double fx = (double)DotProd(v02, v01) / squareLength(v01);
				double fy = (double)DotProd(v02, v01Rot90) / squareLength(v01Rot90);
				triangleLocal[i][j][0] = fx; 
				triangleLocal[i][j][1] = fy; 

				// Sanity check 
				System.out.println("In initializeMesh(). Check validity of fx and fy computation");
				System.out.println("Difference: " + ((v0[0] + fx*v01[0] + fy*v01Rot90[0] - v2[0])*(v0[0] + fx*v01[0] + fy*v01Rot90[0] - v2[0]) +
				(v0[1] + fx*v01[1] + fy*v01Rot90[1] - v2[1])*(v0[1] + fx*v01[1] + fy*v01Rot90[1] - v2[1])));
			}
		}
	}


	public int[] addConstraint(int x, int y) {
		int distToPoint;
		int[] aConstraint; 
		for(int i=0; i<constrainedIdx.size(); i++) {
			aConstraint = deformedVertices[constrainedIdx.get(i)]; 
			distToPoint = (aConstraint[0]-x)*(aConstraint[0]-x) +
				(aConstraint[1]-y)*(aConstraint[1]-y);
			if(distToPoint < 25) {
				return new int[]{0, 0, 0}; 
			}
		}
		for(int i=0; i<numVertices; i++) { 
			distToPoint = (deformedVertices[i][0]-x)*(deformedVertices[i][0]-x) +
				(deformedVertices[i][1]-y)*(deformedVertices[i][1]-y);
			if(distToPoint < 25) {
				constrainedIdx.add(i);
				return new int[]{1, deformedVertices[i][0], deformedVertices[i][1]}; 
			}
		}

		return new int[]{0, 0, 0};
	}

	public void calcStepOneMatrix() { 
		double[][] gMatrix = new double[2*numVertices][2*numVertices];
		// Zero 
		for(int i=0; i<2*numVertices; i++) {
			for(int j=0; j<2*numVertices; j++) {
				gMatrix[i][j] = 0.0;
			}
		}

		// Number of constraints 
		int numConstraints = constrainedIdx.size(); 
		int numFreeVars = numVertices - numConstraints;

		// Figure out vertex ordering 
		vertMap = new HashMap<XYKey, Integer>();

		boolean isConstraint; 
		int nRow = 0; 
		for(int i=0; i<numVertices; i++) {
			isConstraint = false; 
			for(int j=0; j<numConstraints; j++) {
				if(constrainedIdx.get(j) == i) {
					isConstraint = true; 
					break;
				}
			}
			if(!isConstraint) {
				//System.out.println(i);
				//System.out.println(initialVertices[i][0] + ", " + initialVertices[i][1]);
				vertMap.put(new XYKey(initialVertices[i][0],initialVertices[i][1]), nRow);
				nRow++; 
			} 
		}
		for(int i=0; i<numConstraints; i++) {
			//System.out.println(constrainedIdx.get(i));
			//System.out.println(initialVertices[constrainedIdx.get(i)][0] +
			//	", " + initialVertices[constrainedIdx.get(i)][1]);
			vertMap.put(new XYKey(initialVertices[constrainedIdx.get(i)][0], 
				initialVertices[constrainedIdx.get(i)][1]), nRow);
			nRow++; 
		}
		//System.out.println(vertMap.size());
		//System.out.println(initialVertices.length);	

		// Build vector for computing the G matrix 
		vVector = new int[2*numVertices];
		for(int i=0; i<numVertices; i++) {
			isConstraint = false; 
			for(int j=0; j<constrainedIdx.size(); j++) {
				if(i == constrainedIdx.get(j)) {
					isConstraint = true; 
					break; 
				}
			} 
			if(!isConstraint) {
				int rowInVec = vertMap.get(new XYKey(initialVertices[i][0], initialVertices[i][1]));
				vVector[2*rowInVec] = initialVertices[i][0];
				vVector[2*rowInVec + 1] = initialVertices[i][1]; 
			}
		}
		for(int i=0; i<numConstraints; i++) {
			int rowInVec = vertMap.get(
				new XYKey(initialVertices[constrainedIdx.get(i)][0], initialVertices[constrainedIdx.get(i)][1]));
			vVector[2*rowInVec] = initialVertices[constrainedIdx.get(i)][0];
			vVector[2*rowInVec + 1] = initialVertices[constrainedIdx.get(i)][1]; 
		}

		// Build the G Matrix 
		int n0x, n0y, n1x, n1y, n2x, n2y; 
		double x, y; 
		int[] v0 = new int[2], v1 = new int[2], v2 = new int[2];
		int[] v01 = new int[2], v01Rot90 = new int[2];
		int[] v02 = new int[2];
		double v0x, v0y, v1x, v1y, v2x, v2y;
		for(int i=0; i<numTriangles; i++) {
			Polygon triangle = triangles[i]; 
			double triangleError = 0; 
			for(int j=0; j<3; ++j) {
				n0x = 2 * vertMap.get(new XYKey(triangle.xpoints[j], triangle.ypoints[j]));
				n0y = n0x + 1; 
				n1x = 2 * vertMap.get(new XYKey(triangle.xpoints[(j+1)%3], triangle.ypoints[(j+1)%3]));
				n1y = n1x + 1; 
				n2x = 2 * vertMap.get(new XYKey(triangle.xpoints[(j+2)%3], triangle.ypoints[(j+2)%3]));
				n2y = n2x + 1; 
				x = triangleLocal[i][j][0];
				y = triangleLocal[i][j][1];

				
				System.out.println("In calcStepOneMatrix(). Check validity of x and y values.");
				v0[0] = vVector[n0x];
				v0[1] = vVector[n0y];
				v1[0] = vVector[n1x];
				v1[1] = vVector[n1y];
				v2[0] = vVector[n2x];
				v2[1] = vVector[n2y];
				v01[0] = v1[0] - v0[0];
				v01[1] = v1[1] - v0[1];
				v01Rot90[0] = v01[1];
				v01Rot90[1] = -v01[0];
				System.out.println("Error for Vertex: " + ((v0[0] + x*v01[0] + y*v01Rot90[0] - v2[0])
					*(v0[0] + x*v01[0] + y*v01Rot90[0] - v2[0]) + 
					(v0[1] + x*v01[1] + y*v01Rot90[1] - v2[1])
					*(v0[1] + x*v01[1] + y*v01Rot90[1] - v2[1])));
				System.out.println(" ");
				
				System.out.println("In calcStepOneMatrix(). Error equation check.");
				v0x = (double)v0[0];
				v0y = (double)v0[1];
				v1x = (double)v1[0];
				v1y = (double)v1[1];
				v2x = (double)v2[0];
				v2y = (double)v2[1];

				/*
				System.out.println("v0[0]: " + v0[0]);
				System.out.println("v0x: " + v0x);
				System.out.println("v0[1]: " + v0[1]);
				System.out.println("v0y: " + v0y);
				System.out.println("v1[0]: " + v1[0]);
				System.out.println("v1x: " + v1x);
				System.out.println("v1[1]: " + v1[1]);
				System.out.println("v1y: " + v1y);
				System.out.println("v2[0]: " + v2[0]);
				System.out.println("v2x: " + v2x);
				System.out.println("v2[1]: " + v2[1]);
				System.out.println("v2y: " + v2y);
				
				System.out.println("v01[0]: " + v01[0]);
				System.out.println("(v1x-v0x): " + (v1x-v0x));
				System.out.println("v01[1]: " + v01[1]);
				System.out.println("(v1y-v0y): " + (v1y-v0y));
				
				System.out.println("v01Rot90[0]: " + v01Rot90[0]);
				System.out.println("(v1y-v0y): " + (v1y-v0y));
				System.out.println("v01Rot90[1]: " + v01Rot90[1]);
				System.out.println("(v0x-v1x): " + (v0x-v1x));
				*/


				
				double vertexError = 0.0; 
				vertexError += (v0x + x*(v1x-v0x) + y*(v1y-v0y) - v2x)*(v0x + x*(v1x-v0x) + y*(v1y-v0y) - v2x);
				vertexError += (v0y + x*(v1y-v0y) + y*(v0x-v1x) - v2y)*(v0y + x*(v1y-v0y) + y*(v0x-v1x) - v2y);
				System.out.println("Vertex Error: " + vertexError);
				
				/*
				vertexError = 0.0;
				vertexError += v0x*v0x*(x*x - 2*x + y*y + 1);
				vertexError += v0x*v1x*(-2*x*x + 2*x - 2*y*y);
				vertexError += v0x*v0y*(2*y);
				vertexError += v0x*v2x*(2*x - 2);
				vertexError += v0x*v2y*(-2*y);
				vertexError += v0y*v0y*(x*x - 2*x + y*y + 1);
				vertexError += v0y*v1x*(-2*y);
				vertexError += v0y*v1y*(-2*x*x + 2*x - 2*y*y);
				vertexError += v0y*v2x*(2*y);
				vertexError += v0y*v2y*(2*x - 2);
				vertexError += v1x*v1x*(x*x + y*y);
				vertexError += v1x*v2x*(-2*x);
				vertexError += v1x*v2y*(2*y);
				vertexError += v1y*v1y*(x*x + y*y);
				vertexError += v1y*v2x*(-2*y);
				vertexError += v1y*v2y*(-2*x);
				vertexError += v2x*v2x;
				vertexError += v2y*v2y;
				
				System.out.println("Vertex Error: " + vertexError);
				System.out.println(" ");
				*/
				
				/*
				double dTest = 
				(1 - 2*x + (x*x) + (y*y))*(vVector[n0x]*vVector[n0x]) +
				(1 - 2*x + (x*x) + (y*y))*(vVector[n0y]*vVector[n0y]) + 
				((x*x) + (y*y))*(vVector[n1x]*vVector[n1x]) + 
				((x*x) + (y*y))*(vVector[n1y]*vVector[n1y]) + 
				(vVector[n2x]*vVector[n2x]) + (vVector[n2y]*vVector[n2y]) + 
				vVector[n1y]*(-2*y*vVector[n2x] - 2*x*vVector[n2y]) + 
				vVector[n0y]*(-2*y*vVector[n1x] + (2*x - 2*(x*x) - 2*(y*y))*vVector[n1y] + 
				2*y*vVector[n2x] + 
				(-2 + 2*x)*vVector[n2y]) + 
				vVector[n0x]*((2*x - 2*(x*x) - 2*(y*y))*vVector[n1x] + 2*y*vVector[n1y] + (-2 + 2*x)*vVector[n2x] - 
				2*y*vVector[n2y]) + vVector[n1x]*(-2*x*vVector[n2x] + 2*y*vVector[n2y]);
				System.out.println("TEST: " + dTest);
				*/

				
				gMatrix[n0x][n0x] += 1 - 2*x + x*x + y*y;
				gMatrix[n0x][n1x] += 2*x - 2*x*x - 2*y*y;		
				gMatrix[n0x][n1y] += 2*y;						
				gMatrix[n0x][n2x] += -2 + 2*x;					
				gMatrix[n0x][n2y] += -2 * y;	
				gMatrix[n0y][n0y] += 1 - 2*x + x*x + y*y;
				gMatrix[n0y][n1x] += -2*y;						
				gMatrix[n0y][n1y] += 2*x - 2*x*x - 2*y*y;		
				gMatrix[n0y][n2x] += 2*y;						
				gMatrix[n0y][n2y] += -2 + 2*x;	
				gMatrix[n1x][n1x] += x*x + y*y;
				gMatrix[n1x][n2x] += -2*x;						
				gMatrix[n1x][n2y] += 2*y;	
				gMatrix[n1y][n1y] += x*x + y*y;
				gMatrix[n1y][n2x] += -2*y;						
				gMatrix[n1y][n2y] += -2*x;	
				gMatrix[n2x][n2x] += 1;
				gMatrix[n2y][n2y] += 1;

				triangleError += (1 - 2*x + x*x + y*y)  * vVector[n0x] * vVector[n0x];
				triangleError += (2*x - 2*x*x - 2*y*y)  * vVector[n0x] * vVector[n1x];
				triangleError += (2*y)                  * vVector[n0x] * vVector[n1y];
				triangleError += (-2 + 2*x )            * vVector[n0x] * vVector[n2x];
				triangleError += (-2 * y)               * vVector[n0x] * vVector[n2y];
				triangleError += (1 - 2*x + x*x + y*y)   * vVector[n0y] * vVector[n0y];
				triangleError += (-2*y)                  * vVector[n0y] * vVector[n1x];
				triangleError += (2*x - 2*x*x - 2*y*y)   * vVector[n0y] * vVector[n1y];
				triangleError += (2*y)                   * vVector[n0y] * vVector[n2x];
				triangleError += (-2 + 2*x)              * vVector[n0y] * vVector[n2y];
				triangleError += (x*x + y*y)            * vVector[n1x] * vVector[n1x];
				triangleError += (-2*x)                 * vVector[n1x] * vVector[n2x];
				triangleError += (2*y)                  * vVector[n1x] * vVector[n2y];
				triangleError += (x*x + y*y)            * vVector[n1y] * vVector[n1y];
				triangleError += (-2*y)                 * vVector[n1y] * vVector[n2x];
				triangleError += (-2*x)                 * vVector[n1y] * vVector[n2y];
				triangleError += vVector[n2x] * vVector[n2x]  +  vVector[n2y] * vVector[n2y] ;
				
			}
			System.out.println("Triangle Error: " + triangleError);
		}
		// Sanity Check 
		/*
		Matrix sanityG = new Matrix(gMatrix);
		double sanityV[][] = new double[vVector.length][1];
		for(int i=0; i<vVector.length; i++) {
			sanityV[i][0] = vVector[i];
		}
		Matrix sanityVMatrix = new Matrix(sanityV);
		Matrix error = sanityVMatrix.transpose().times(sanityG.times(sanityVMatrix));
		System.out.println(" ");
		System.out.println("Total Error from Matrix Mult: " + error.getArrayCopy()[0][0]);
		System.out.println(" ");
		*/

		System.out.println("Num Free Var: " + numFreeVars);
		System.out.println("Num Constraints: " + numConstraints);
		Matrix g = new Matrix(gMatrix);
		Matrix g00 = g.getMatrix(0, 2*numFreeVars-1, 0, 2*numFreeVars-1);
		Matrix g01 = g.getMatrix(0, 2*numFreeVars-1, 2*numFreeVars, 2*numVertices-1);
		Matrix g10 = g.getMatrix(2*numFreeVars, 2*numVertices-1, 0, 2*numFreeVars-1);
		// Compute GPrime 
		Matrix gPrime = g00.plus(g00.transpose());
		Matrix b = g01.plus(g10.transpose());
		Matrix gPrimeInv = gPrime.inverse();
		gPrimeInvB = gPrimeInv.times(b);
		gPrimeInvB = gPrimeInvB.times(-1);


		// Sanity check 
		/*
		System.out.println("Final Matrix Dimensions: " + finalMatrix.getRowDimension() + 
			", " + finalMatrix.getColumnDimension());
		// Constraints 
		double[][] qVec = new double[2*numConstraints][1];
		for(int i=0; i<2*numConstraints; i++) {
			qVec[i][0] = vVector[2*numFreeVars + i];
		}
		Matrix qVecMat = new Matrix(qVec);
		Matrix results = finalMatrix.times(qVecMat);
		System.out.println("Results dimensions: " + 
			results.getRowDimension() + 
			", " + results.getColumnDimension());
		for(int i=0; i<2*numFreeVars; i+=2) {
			System.out.println("Original: " + 
				vVector[i] + ", " + 
				vVector[i+1]); 
			System.out.println("Result: " + 
				results.get(i,0)+ ", " + 
				results.get(i+1,0)); 
		}
		*/
	}
	

	
	public int[] setConstraintActive(int x, int y) {
		int distToPoint; 
		int numConstraints = constrainedIdx.size();
		for(int i=0; i<numConstraints; i++) {
			int aConstraintX = deformedVertices[constrainedIdx.get(i)][0];
			int aConstraintY = deformedVertices[constrainedIdx.get(i)][1];
			distToPoint = (aConstraintX - x)*(aConstraintX - x) + (aConstraintY - y)*(aConstraintY - y); 
			if(distToPoint < 25) { 
				return new int[]{1, constrainedIdx.get(i)};
			} 
		}
		return new int[]{0, 0};
	}
	
	
	// Calculating new coordinates 
	public void stepOne(int idx, int newX, int newY) {
		System.out.println("Constraint Old: " + deformedVertices[idx][0] + ", " + 
			deformedVertices[idx][1]);
		System.out.println("Constraint New: " + newX + ", " + newY);
		// Update deformed vertex 
		deformedVertices[idx][0] = newX;
		deformedVertices[idx][1] = newY; 
		// Compile qVector 
		int numConstraints = constrainedIdx.size();
		int numFreeVars = numVertices - numConstraints;
		double[][] qVec = new double[2*numConstraints][1];
		for(int i=0; i<numConstraints; i++) {
			qVec[2*i][0] = deformedVertices[constrainedIdx.get(i)][0];
			qVec[2*i + 1][0] = deformedVertices[constrainedIdx.get(i)][1];
		}
		// Sanity Check 
		/* 
		for(int i=0; i<2*numConstraints; i++) {
			System.out.println(qVec[i][0] + ", " + vVector[2*numFreeVars + i]);
		}
		*/ 
		// Build map to initial vertices 
		HashMap<XYKey, Integer> mapToInitial = new HashMap<XYKey, Integer>(); 
		for(int i=0; i < numVertices; i++) { 
			mapToInitial.put(new XYKey(initialVertices[i][0], initialVertices[i][1]), i); 
		}
		Matrix qVecMat = new Matrix(qVec);
		Matrix results = gPrimeInvB.times(qVecMat);
		int origX, origY, origIdx; 
		for(int i=0; i<2*numFreeVars; i+=2) {
			origX = vVector[i]; 
			origY = vVector[i+1]; 
			origIdx = mapToInitial.get(new XYKey(origX, origY)); 
			if(origIdx != -1) {
 				deformedVertices[origIdx][0] = (int)results.get(i, 0); 
				deformedVertices[origIdx][1] = (int)results.get(i+1, 0);
			}
		}
		// Sanity Check 
		for(int i=0; i<numVertices; i++) {
			System.out.println("Original: (" + initialVertices[i][0] +
				", " + initialVertices[i][1] + "), " +
				"Final: (" + deformedVertices[i][0] + ", " + 
				deformedVertices[i][1] + ")");
		}
		// Update Polygon 
		updatedTri = new Polygon[numTriangles];
		int[] newTriX = new int[3];
		int[] newTriY = new int[3];
		int currX, currY, currIdx; 
		for(int i=0; i<numTriangles; i++) {
			for(int j=0; j<3; j++){
				currX = triangles[i].xpoints[j];
				currY = triangles[i].ypoints[j]; 
				currIdx = mapToInitial.get(new XYKey(currX, currY)); 
				newTriX[j] = deformedVertices[currIdx][0];
				newTriY[j] = deformedVertices[currIdx][1];
			}
			updatedTri[i] = new Polygon(newTriX, newTriY, 3);
		}
		int[] newPolyX = new int[numVertices]; 
		int[] newPolyY = new int[numVertices];
		for(int i=0; i<numVertices; i++) { 
			currX = xpoints[i];
			currY = ypoints[i]; 
			currIdx = mapToInitial.get(new XYKey(currX, currY));
			newPolyX[i] = deformedVertices[currIdx][0];
			newPolyY[i] = deformedVertices[currIdx][1];
		}
		updatedPoly = new Polygon(newPolyX, newPolyY, numVertices);
	}
}