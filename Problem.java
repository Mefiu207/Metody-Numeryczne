import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


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

    ArrayList<Double[][]> u;

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

        // Initialize ArrayList
        u = new ArrayList<Double[][]>();           
    }

    public readPsiFromFile(String filename){
        
        try{    
            File myFile = new File(filename);
            Scanner myScanner = new Scanner(myFile);

            while()
        } catch(FileNotFoundException e) { e.printStackTrace(); }

    }


}