//class CCF1D : public CorrFunction {
//
//    //convencao 
//    //target[0][i] = par 0 0
//    //target[1][i] = par 0 1
//    //target[2][i] = par 1 0
//    //target[3][i] = par 1 1
//
//public:
//    vector<vector<int>> aux = { {}, {}, {}, {} };
//    vector<vector<int>> curr = { {}, {}, {}, {} };
//    vector<vector<double>> target = { {}, {}, {}, {} };
//    Microstructure& mic;
//    int axis, rmax;
//    double currEnergy, auxEnergy;
//
//
//
//    CCF1D(Microstructure& mic_, string targetPath, int axis_) : mic(mic_) {
//        axis = axis_;
//
//        if (axis == 0) {
//            rmax = mic.nx / 2;
//        }
//        else if (axis == 1) {
//            rmax = mic.ny / 2;
//        }
//        else if (axis == 2) {
//            rmax = mic.nz / 2;
//        }
//
//        firstCalc();
//        readFile(targetPath);
//        currEnergy = getStartEnergy();
//
//    }
//
//    void firstCalc() {
//        curr[0].assign(rmax, 0);
//        curr[1].assign(rmax, 0);
//        curr[2].assign(rmax, 0);
//        curr[3].assign(rmax, 0);
//
//        for (int i = 0; i < mic.n; i++) {
//            mic.setPoint0(i);
//            mic.setPoint1(i);
//            for (int r = 0; r < rmax; r++) {
//                if (mic.phasePoint0() == 0 and mic.phasePoint1() == 0) curr[0][r]++;
//                if (mic.phasePoint0() == 0 and mic.phasePoint1() == 1) curr[1][r]++;
//                if (mic.phasePoint0() == 1 and mic.phasePoint1() == 0) curr[2][r]++;
//                if (mic.phasePoint0() == 1 and mic.phasePoint1() == 1) curr[3][r]++;
//                mic.advancePoint1(axis);
//            }
//
//        }
//    }
//
//    double getStartEnergy() {
//        double sum = 0;
//        for (int r = 0; r < rmax; r++) {
//            for (int i = 0; i < 4; i++) {
//                sum += pow(target[i][r] - (double)curr[i][r] / mic.n, 2);
//            }
//        }
//        return sum;
//    }
//
//    double getAuxEnergy() override {
//        double sum = 0;
//        for (int r = 0; r < rmax; r++) {
//            for (int i = 0; i < 4; i++) {
//                sum += pow(target[i][r] - (double)aux[i][r] / mic.n, 2);
//            }
//        }
//        return sum;
//    }
//
//
//    double energyDiff() override {
//        return auxEnergy - currEnergy;
//    }
//
//    void swap() override {
//        curr = aux;
//        currEnergy = auxEnergy;
//    }
//
//
//    void readFile(string path) {
//
//        ifstream file(path);
//        string trash;
//        string value00;
//        string value01;
//        string value10;
//        string value11;
//
//        string line;
//
//        target[0].assign(rmax, 0);
//        target[1].assign(rmax, 0);
//        target[2].assign(rmax, 0);
//        target[3].assign(rmax, 0);
//
//        for (int i = 0; i < rmax; i++) {
//            file >> trash >> value00 >> value01 >> value10 >> value11;
//            target[0][i] = stof(value00);
//            target[1][i] = stof(value01);
//            target[2][i] = stof(value10);
//            target[3][i] = stof(value11);
//            getline(file, line);
//        }
//
//        file.close();
//    }
//
//    void update(int n0_, int n1_) override {
//        //p0_ é sempre da fase 1
//        //p1_ é sempre da fase 2
//
//        aux = curr;
//
//        mic.setPoint0(n0_);
//        mic.setPoint1(n0_);
//
//        //p0
//        for (int r = 1; r < rmax; r++) {
//            mic.backPoint0(axis);
//            mic.advancePoint1(axis);
//
//            if (mic.phasePoint0() == 0) {
//                aux[0][r]--;
//                aux[1][r]++;
//            }
//            else if (mic.phasePoint0() == 1) {
//                aux[2][r]--;
//                aux[3][r]++;
//            }
//
//            if (mic.phasePoint1() == 0) {
//                aux[0][r]--;
//                aux[2][r]++;
//            }
//            else if (mic.phasePoint1() == 1) {
//                aux[1][r]--;
//                aux[3][r]++;
//            }
//        }
//
//        //p1
//        mic.setPoint0(n1_);
//        mic.setPoint1(n1_);
//
//        for (int r = 1; r < rmax; r++) {
//            mic.backPoint0(axis);
//            mic.advancePoint1(axis);
//
//            if (mic.phasePoint0() == 0) {
//                aux[0][r]++;
//                aux[1][r]--;
//
//            }
//            else if (mic.phasePoint0() == 1) {
//                aux[2][r]++;
//                aux[3][r]--;
//            }
//
//            if (mic.phasePoint1() == 0) {
//                aux[0][r]++;
//                aux[2][r]--;
//            }
//            else if (mic.phasePoint1() == 1) {
//                aux[1][r]++;
//                aux[3][r]--;
//            }
//        }
//
//        auxEnergy = getAuxEnergy();
//    }
//
//    double getCurrEnergy() override {
//        return currEnergy;
//    }
//
//
//};
//
//class TwoPoint1D : public CorrFunction {
//
//public:
//    vector<int> aux;
//    vector<int> curr;
//    vector<double> target;
//    Microstructure& mic;
//    int axis, rmax, phase;
//    double currEnergy;
//    double auxEnergy;
//
//
//
//    TwoPoint1D(Microstructure& mic_, string targetPath, int axis_, int phase_) : mic(mic_) {
//        axis = axis_;
//        phase = phase_;
//
//        if (axis == 0) {
//            rmax = mic.nx / 2;
//        }
//        else if (axis == 1) {
//            rmax = mic.ny / 2;
//        }
//        else if (axis == 2) {
//            rmax = mic.nz / 2;
//        }
//
//        firstCalc();
//        readFile(targetPath);
//        currEnergy = getStartEnergy();
//
//    }
//
//
//    void readFile(string path) {
//
//        ifstream file(path);
//        string value;
//
//        string line;
//
//        target.assign(rmax, 0);
//
//        for (int i = 0; i < rmax; i++) {
//            file >> value >> value;
//            target[i] = stof(value);
//            getline(file, line);
//        }
//
//        file.close();
//    }
//
//    void firstCalc() {
//        curr.assign(rmax, 0);
//
//        for (int i = 0; i < mic.n; i++) {
//            mic.setPoint0(i);
//            mic.setPoint1(i);
//            for (int r = 0; r < rmax; r++) {
//                if (mic.phasePoint0() == phase and mic.phasePoint1() == phase) curr[r] = curr[r] + 1;
//                mic.advancePoint1(axis);
//            }
//
//        }
//    }
//
//    void update(int n0_, int n1_) override {
//        //p0_ é sempre da fase 1
//        //p1_ é sempre da fase 2
//
//        aux = curr;
//
//        mic.setPoint0(n0_);
//        mic.setPoint1(n0_);
//
//        //p0
//        for (int r = 1; r < rmax; r++) {
//            mic.backPoint0(axis);
//            mic.advancePoint1(axis);
//
//            if (mic.phasePoint1() == phase) aux[r]++;
//            if (mic.phasePoint0() == phase) aux[r]++;
//        }
//
//        //p1
//        mic.setPoint0(n1_);
//        mic.setPoint1(n1_);
//
//        for (int r = 1; r < rmax; r++) {
//            mic.backPoint0(axis);
//            mic.advancePoint1(axis);
//
//            if (mic.phasePoint1() == phase) aux[r]--;
//            if (mic.phasePoint0() == phase) aux[r]--;
//
//        }
//
//        auxEnergy = getAuxEnergy();
//    }
//
//    double energyDiff() override {
//        return auxEnergy - currEnergy;
//    }
//
//    void swap() override {
//        curr = aux;
//        currEnergy = auxEnergy;
//    }
//
//    double getAuxEnergy() override {
//        double sum = 0;
//        for (int i = 0; i < rmax; i++) {
//            sum += pow(target[i] - (double)aux[i] / mic.n, 2);
//        }
//        return sum;
//    }
//
//    double getCurrEnergy() override {
//        return currEnergy;
//    }
//
//    double getStartEnergy() {
//        double sum = 0;
//        for (int i = 0; i < rmax; i++) {
//            sum += pow(target[i] - (double)curr[i] / mic.n, 2);
//        }
//        return sum;
//    }
//
//};

//class Point {
//
//public:
//    int x = 0;
//    int y = 0;
//    int z = 0;
//    Microstructure& mic;
//
//    Point(Microstructure& mic_, int x_, int y_, int z_) : mic(mic_) {
//        setPoint(x_, y_, z_);
//    }
//
//    Point(Microstructure& mic_, int n) : mic(mic_) {
//        setPoint(n);
//    }
//
//    Point(Microstructure& mic_) : mic(mic_) {
//    }
//
//
//    void correctX() {
//        if (x >= mic.nx) { x = x - mic.nx; }
//        else if (x < 0) { x = x + mic.nx; };
//    }
//
//    void correctY() {
//        if (y >= mic.ny) { y = y - mic.ny; }
//        else if (y < 0) { y = y + mic.ny; }
//    }
//
//    void correctZ() {
//        if (z >= mic.nz) { z = z - mic.nz; }
//        else if (z < 0) { z = z + mic.nz; }
//    }
//
//    void advance(int axis, int dist) {
//        switch (axis) {
//        case(0):
//            x = x + dist;
//            correctX();
//            break;
//        case(1):
//            y = y + dist;
//            correctY();
//            break;
//        case(2):
//            z = z + dist;
//            correctZ();
//            break;
//        }
//    }
//
//    void next(int axis) {
//        advance(axis, 1);
//    }
//
//    void prev(int axis) {
//        advance(axis, -1);
//    }
//
//    int setPoint(int xx, int yy, int zz) {
//        x = xx;
//        y = yy;
//        z = zz;
//    }
//
//    void setPoint(int n) {
//        z = n / (mic.nx * mic.ny);
//        y = n % (mic.nx * mic.ny);
//        x = y % mic.nx;
//        y = y / mic.nx;
//    }
//
//    int ind() {
//        return mic.nx * mic.ny * z + mic.nx * y + x;
//    }
//    int getPhase() {
//        return mic.arr[ind()];
//    }
//
//};