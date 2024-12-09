import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.io.FileWriter;
import java.io.IOException;

public class Problem{
    
    Double dt;
    Double delta;

    Integer nx;
    Integer ny;

    Double[] x;
    Double[] y;

    Integer i1;
    Integer i2;
    Integer j1;

    Double sigma;

    Double xa;
    Double ya;

    Double psi[][];

    Double vx[][];
    Double vy[][];

    Double u0[][];
    Double u1[][];

    public Problem(){
        delta = 0.01;

        nx = 400;
        ny = 90;

        i1 = 200;
        i2 = 210;
        j1 = 50;

        sigma = 10 * delta;

        xa = 0.45;
        ya = 0.45;

        // Initialize arrays
        x = new Double[nx + 1];
        y = new Double[ny + 1];

        psi = new Double[ny + 1][nx + 1];
        vx = new Double[ny + 1][nx + 1];
        vy = new Double[ny + 1][nx + 1];

        u0 = new Double[ny + 1][nx + 1];
        u1 = new Double[ny + 1][nx + 1];

    }

    // Reades data from a given file to psi array
    public void readPsiFromFile(String filename){
        
        try{    
            File myFile = new File(filename);
            Scanner myScanner = new Scanner(myFile);

            while(myScanner.hasNextLine()){

                String line = myScanner.nextLine().trim();

                String[] parts = line.split("\\s+");

                Integer xPoint = Integer.parseInt(parts[0]);
                Integer yPoint = Integer.parseInt(parts[1]);
                
                Double value = Double.parseDouble(parts[2]);

                psi[yPoint][xPoint] = value;

                myScanner.close();
            }

        } catch(FileNotFoundException e) { e.printStackTrace(); }

    }

    // Calculates values for Vx and Vy arrays
    public void calcV(){
        
        // Inner field
        for(int j = 1; j < ny; j++){
            for(int i = 1; i < nx; i++){
                
                vx[j][i] = (psi[j + 1][i] - psi[j - 1][i]) / (2 * delta);
                vy[j][i] = -(psi[j][i + 1] - psi[j][i - 1]) / (2 * delta);

            }
        }

        // Edges inside the box
        for(int j = 0; j <= j1; j++){
            for(int i = i1; i <= i2; i++){

                vx[j][i] = 0.0;
                vy[j][i] = 0.0;

            }
        }

        // Upper and lower edge
        for(int i = 1; i < nx; i++){

            vx[0][i] = 0.0;
            vy[ny][i] = 0.0;

        }

        // Left and right edge
        for(int j = 0; j <= ny; j++){

            vx[j][0] = vx[j][1];
            vx[j][nx] = vx[j][nx - 1];

        }

        setDeltaT();

    }    

   // Finds Vmax and sets value of dt (delta T) (calulated in calV() method)
    private void setDeltaT(){
        Double currentHighest = 0.0;

        Double presentV = 0.0;

        for(int j = 1; j < ny; j++){
            for(int i = 1; i < nx; i++){

                presentV = Math.sqrt(Math.pow(vx[j][i], 2) + Math.pow(vy[j][i], 2)); 
                if(presentV > currentHighest) currentHighest = presentV;

            }
        }

        this.dt = delta/(4*currentHighest);
    }

    // Saves values of vx and vy to given file
    public void saveVToFile(String filename){
        
        File myFile = new File(filename);
        myFile.delete();

        try {myFile.createNewFile();}
        catch(IOException e) {e.printStackTrace();}

        try{
            FileWriter writeFile = new FileWriter(filename);

            for(int j = 0; j <= ny; j++){
                for(int i = 0; i <= nx; i++){
                    writeFile.write(i + " " + j + " " + vx[j][i] + " " + vy[j][i] + "\n");
                }
                writeFile.write("\n");
            }
        }catch(IOException e) { e.printStackTrace(); }   
    } 

    // Calculates the first u0 value 
    private void setU0(){

        Double firstPart = 1/(2*Math.PI*Math.pow(sigma, 2));

        for(int j = 0; j <= ny; j++){
            for(int i = 0; i <= nx; i++){

                u0[j][i] = firstPart * Math.exp(-(Math.pow(x[i] - xa, 2) + Math.pow(y[j] - ya, 2)) / (2* Math.pow(sigma, 2))); 

            }
        }

    }

    public void picardIter(){
        setU0();

    }





    public static void main(String[] args){
        
        Problem a1 = new Problem();
        
        a1.readPsiFromFile("psi.dat");
        a1.calcV();
        a1.saveVToFile("v.dat");
    
    }
}