//class CCF : public CorrFunction {
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
//    int rmax, numAxis;
//    double currEnergy, auxEnergy;
//    bool isReference;
//
//    CCF(Microstructure& mic_, int numAxis_) : mic(mic_) { //lê mic e assume que é referência (sem target)
//        isReference = true;
//
//        numAxis = numAxis_;
//        if (numAxis == 1) {
//            rmax = mic.nx / 2;
//        }
//        else if (numAxis == 2) {
//            rmax = min(mic.nx / 2, mic.ny / 2);
//        }
//        else if (numAxis == 3) {
//            rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
//        }
//
//        firstCalc();
//
//        target[0].assign(rmax, 0);
//        target[1].assign(rmax, 0);
//        target[2].assign(rmax, 0);
//        target[3].assign(rmax, 0);
//
//        for (int r = 0; r < rmax; r++) {
//            for (int i = 0; i < 4; i++) {
//                target[i][r] = double(curr[i][r]) / (numAxis * mic.n);
//            }
//        }
//
//
//    }
//
//
//    CCF(Microstructure& mic_, vector<vector<double>> target_, int numAxis_) : mic(mic_) { //aceita vetor como entrada
//        isReference = false;
//
//        numAxis = numAxis_;
//        if (numAxis == 1) {
//            rmax = mic.nx / 2;
//        }
//        else if (numAxis == 2) {
//            rmax = min(mic.nx / 2, mic.ny / 2);
//        }
//        else if (numAxis == 3) {
//            rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
//        }
//
//        firstCalc();
//        target = target_;
//        currEnergy = getStartEnergy();
//
//    }
//
//    CCF(Microstructure& mic_, string targetPath, int numAxis_) : mic(mic_) { // lê função do arquivo
//        isReference = false;
//
//        numAxis = numAxis_;
//        if (numAxis == 1) {
//            rmax = mic.nx / 2;
//        }
//        else if (numAxis == 2) {
//            rmax = min(mic.nx / 2, mic.ny / 2);
//        }
//        else if (numAxis == 3) {
//            rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
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
//        for (int axis = 0; axis < numAxis; axis++) {
//            for (int i = 0; i < mic.n; i++) {
//                mic.setPoint0(i);
//                mic.setPoint1(i);
//                for (int r = 0; r < rmax; r++) {
//                    if (mic.phasePoint0() == 0 and mic.phasePoint1() == 0) curr[0][r]++;
//                    if (mic.phasePoint0() == 0 and mic.phasePoint1() == 1) curr[1][r]++;
//                    if (mic.phasePoint0() == 1 and mic.phasePoint1() == 0) curr[2][r]++;
//                    if (mic.phasePoint0() == 1 and mic.phasePoint1() == 1) curr[3][r]++;
//                    mic.advancePoint1(axis);
//                }
//            }
//        }
//    }
//
//    double getStartEnergy() {
//        double sum = 0;
//        for (int r = 0; r < rmax; r++) {
//            for (int i = 0; i < 4; i++) {
//                sum += pow(target[i][r] - double(curr[i][r]) / (numAxis * mic.n), 2);
//            }
//        }
//        return sum;
//    }
//
//    double getAuxEnergy() override {
//        double sum = 0;
//        for (int r = 0; r < rmax; r++) {
//            for (int i = 0; i < 4; i++) {
//                sum += pow(target[i][r] - double(aux[i][r]) / (numAxis * mic.n), 2);
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
//        //troca
//        mic.arr[n0_] = 1;
//        mic.arr[n1_] = 0;
//
//        for (int axis = 0; axis < numAxis; axis++) {
//
//            mic.setPoint0(n0_);
//            mic.setPoint1(n0_);
//
//            for (int r = 1; r < rmax; r++) {
//                mic.backPoint0(axis);
//                mic.advancePoint1(axis);
//                if (mic.phasePoint0() == 0) { //00 -> 01
//                    aux[0][r]--;
//                    aux[1][r]++;
//                }
//                else { //10 -> 11
//                    aux[2][r]--;
//                    aux[3][r]++;
//                }
//
//                if (mic.phasePoint1() == 0) { //00 -> 10
//                    aux[0][r]--;
//                    aux[2][r]++;
//                }
//                else { //01 -> 11
//                    aux[1][r]--;
//                    aux[3][r]++;
//                }
//            }
//
//            mic.setPoint0(n1_);
//            mic.setPoint1(n1_);
//
//            for (int r = 1; r < rmax; r++) {
//                mic.backPoint0(axis);
//                mic.advancePoint1(axis);
//                if (mic.phasePoint0() == 0) { //01 -> 00
//                    aux[1][r]--;
//                    aux[0][r]++;
//                }
//                else { //11 -> 10
//                    aux[3][r]--;
//                    aux[2][r]++;
//                }
//
//                if (mic.phasePoint1() == 0) { //10 -> 00
//                    aux[2][r]--;
//                    aux[0][r]++;
//                }
//                else { //11 -> 01
//                    aux[3][r]--;
//                    aux[1][r]++;
//                }
//            }
//
//        }
//
//        //destroca
//        mic.arr[n0_] = 0;
//        mic.arr[n1_] = 1;
//
//        auxEnergy = getAuxEnergy();
//    }
//
//    double getCurrEnergy() override {
//        return currEnergy;
//    }
//
//    void writeToFile(string path_) {
//        ofstream file;
//        file.open(path_);
//        double den = mic.n * numAxis;
//        int pos = 0;
//
//        for (int r = 0; r < rmax; r++)
//        {
//            file << r << "\t" << curr[0][r] / den << "\t" << curr[1][r] / den << "\t" << curr[2][r] / den << "\t" << curr[3][r] / den << "\n";
//            pos++;
//        }
//        file.close();
//    }
//
//    void writeToFileNonPeriodic(string path_) {
//        ofstream file;
//        vector<vector<double>> curr = evalNonPeriodic();
//        file.open(path_);
//        int pos = 0;
//
//        for (int r = 0; r < rmax; r++)
//        {
//            file << r << "\t" << curr[0][r] << "\t" << curr[1][r] << "\t" << curr[2][r] << "\t" << curr[3][r] << "\n";
//            pos++;
//        }
//        file.close();
//    }
//
//    vector<vector<double>> evalNonPeriodic() {
//
//        vector<vector<double>> nonPeriodic = { {},{},{},{} };
//        vector<double> den;
//        for (int i = 0; i < 4; i++) nonPeriodic[i].assign(rmax, 0);
//        den.assign(rmax, 0);
//
//        int diff = 0;
//        vector<int> sizes = { mic.nx, mic.ny, mic.nz };
//
//        for (int axis = 0; axis < numAxis; axis++) {
//
//            for (int i = 0; i < mic.n; i++) {
//                mic.setPoint0(i);
//                mic.setPoint1(i);
//                diff = 0;
//                while (diff < rmax) {
//
//                    if (mic.phasePoint0() == 0 and mic.phasePoint1() == 0) nonPeriodic[0][diff]++;
//                    else if (mic.phasePoint0() == 0 and mic.phasePoint1() == 1) nonPeriodic[1][diff]++;
//                    else if (mic.phasePoint0() == 1 and mic.phasePoint1() == 0) nonPeriodic[2][diff]++;
//                    else if (mic.phasePoint0() == 1 and mic.phasePoint1() == 1) nonPeriodic[3][diff]++;
//
//                    mic.advancePoint1(axis);
//                    den[diff]++;
//                    diff++;
//
//                    if (mic.coordPoint1(axis) == mic.axisSize(axis) - 1) break;
//                }
//            }
//        }
//
//        for (int r = 0; r < rmax; r++) {
//            for (int i = 0; i < 4; i++) {
//                nonPeriodic[i][r] = nonPeriodic[i][r] / den[r];
//            }
//        }
//
//        return nonPeriodic;
//
//    }
//
//};
//
//
//class LinealPath : public CorrFunction {
//
//public:
//    vector<int> aux;
//    vector<int> curr;
//    vector<double> target;
//
//    Microstructure& mic;
//    int rmax, numAxis;
//    double currEnergy, auxEnergy;
//    bool isReference = false;
//
//
//    LinealPath(Microstructure& mic_, int numAxis_) : mic(mic_) { //assume como referência
//        isReference = true;
//
//        numAxis = numAxis_;
//
//        if (numAxis == 1) {
//            rmax = mic.nx / 2;
//        }
//        else if (numAxis == 2) {
//            rmax = min(mic.nx / 2, mic.ny / 2);
//        }
//        else if (numAxis == 3) {
//            rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
//        }
//
//        firstCalc();
//
//        target.assign(rmax, 0);
//
//
//
//        for (int r = 0; r < rmax; r++) {
//            target[r] = double(curr[r]) / (numAxis * mic.n);
//        }
//
//
//        currEnergy = getStartEnergy();
//
//    }
//
//    LinealPath(Microstructure& mic_, vector<double> target_, int numAxis_) : mic(mic_) { //lê do vetor
//        numAxis = numAxis_;
//
//        if (numAxis == 1) {
//            rmax = mic.nx / 2;
//        }
//        else if (numAxis == 2) {
//            rmax = min(mic.nx / 2, mic.ny / 2);
//        }
//        else if (numAxis == 3) {
//            rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
//        }
//
//
//        firstCalc();
//
//        target = target_;
//        currEnergy = getStartEnergy();
//
//    }
//
//
//    LinealPath(Microstructure& mic_, string targetPath, int numAxis_) : mic(mic_) { //lê do arquivo
//        numAxis = numAxis_;
//
//        if (numAxis == 1) {
//            rmax = mic.nx / 2;
//        }
//        else if (numAxis == 2) {
//            rmax = min(mic.nx / 2, mic.ny / 2);
//        }
//        else if (numAxis == 3) {
//            rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
//        }
//
//
//        firstCalc();
//
//
//        readFile(targetPath);
//        currEnergy = getStartEnergy();
//
//    }
//
//    void firstCalc() {
//
//        curr.assign(rmax, 0);
//        int lmax = 0;
//
//        for (int i = 0; i < mic.n; i++) {
//            mic.setPoint0(i);
//            if (mic.phasePoint0() == 1) {
//                for (int axis = 0; axis < numAxis; axis++) {
//                    mic.setPoint0(i);
//                    lmax = 0;
//                    mic.advancePoint0(axis);
//
//                    while (mic.phasePoint0() == 1) {
//                        lmax++;
//                        mic.advancePoint0(axis);
//                    }
//
//                    for (int r = 0; r <= lmax; r++) curr[r]++;
//
//                }
//            }
//        }
//    }
//
//    double getStartEnergy() {
//        double sum = 0;
//        for (int r = 0; r < rmax; r++) {
//            sum += pow(target[r] - double(curr[r]) / (numAxis * mic.n), 2);
//        }
//        return sum;
//    }
//
//    double getAuxEnergy() override {
//        double sum = 0;
//        for (int r = 0; r < rmax; r++) {
//            sum += pow(target[r] - double(aux[r]) / (numAxis * mic.n), 2);
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
//        string value;
//
//        string line;
//
//        target.assign(rmax, 0);
//
//        for (int r = 0; r < rmax; r++) {
//            file >> value >> value;
//            target[r] = stod(value);
//            getline(file, line);
//        }
//
//        file.close();
//    }
//
//    void update(int n0_, int n1_) override {
//        //p0_ é sempre da fase 0
//        //p1_ é sempre da fase 1 (fase de interesse)
//
//        aux = curr;
//
//
//        //atualização do tipo 111011 -> 1111111
//
//        int left = 0;
//        int right = 0;
//        int total = 0;
//
//        for (int axis = 0; axis < numAxis; axis++) {
//
//            mic.setPoint0(n0_);
//            mic.setPoint1(n0_);
//            left = 0;
//            right = 0;
//
//            mic.advancePoint1(axis);
//            while (mic.phasePoint1() == 1) {
//                right++;
//                mic.advancePoint1(axis);
//            }
//
//            mic.backPoint0(axis);
//            while (mic.phasePoint0() == 1) {
//                left++;
//                mic.backPoint0(axis);
//            }
//
//            total = left + right + 1; //remover corda de tamanho total
//
//            for (int i = 0; i < left; i++) aux[i] -= (left - i);
//            for (int i = 0; i < right; i++) aux[i] -= (right - i);
//            for (int i = 0; i < total; i++) aux[i] += (total - i);
//
//            //atualização do tipo 111111 -> 1110111
//
//            mic.setPoint0(n1_);
//            mic.setPoint1(n1_);
//            left = 0;
//            right = 0;
//
//            mic.advancePoint1(axis);
//            while (mic.phasePoint1() == 1) {
//                right++;
//                mic.advancePoint1(axis);
//            }
//
//            mic.backPoint0(axis);
//            while (mic.phasePoint0() == 1) {
//                left++;
//                mic.backPoint0(axis);
//            }
//
//            total = left + right + 1; //remover corda de tamanho total
//
//            for (int i = 0; i < total; i++) aux[i] -= (total - i);
//            for (int i = 0; i < left; i++) aux[i] += (left - i);
//            for (int i = 0; i < right; i++) aux[i] += (right - i);
//
//        }
//
//        auxEnergy = getAuxEnergy();
//
//    }
//
//    double getCurrEnergy() override {
//        return currEnergy;
//    }
//
//    void writeToFile(string path_) {
//        ofstream file;
//        file.open(path_);
//        int pos = 0;
//
//        for (int r = 0; r < rmax; r++)
//        {
//            file << r << "\t" << double(curr[r]) / (numAxis * mic.n) << "\n";
//            pos++;
//        }
//        file.close();
//    }
//
//    void writeToFileNonPeriodic(string path_) {
//        ofstream file;
//        vector<double> vec = evalNonPeriodic();
//        file.open(path_);
//        int pos = 0;
//
//        for (int r = 0; r < rmax; r++)
//        {
//            file << r << "\t" << vec[r] << "\n";
//            pos++;
//        }
//        file.close();
//    }
//
//    vector<double> evalNonPeriodic() {
//
//        vector<double> nonPeriodic;
//        vector<int> den;
//        nonPeriodic.assign(rmax, 0);
//        den.assign(rmax, 0);
//
//        for (int r = 0; r < rmax; r++) {
//            if (numAxis == 1) den[r] = mic.n - mic.nx * r;
//            else if (numAxis == 2) den[r] = 2 * mic.n - (mic.nx + mic.ny) * r;
//            else if (numAxis == 3) den[r] = 3 * mic.n - r * (mic.nx * mic.ny + mic.nx * mic.nz + mic.ny * mic.nz);
//        }
//
//
//        int lmax = 0;
//
//        for (int i = 0; i < mic.n; i++) {
//            mic.setPoint0(i);
//            if (mic.phasePoint0() == 1) {
//                for (int axis = 0; axis < numAxis; axis++) {
//                    mic.setPoint0(i);
//                    lmax = 0;
//                    mic.advancePoint0(axis);
//
//                    while (mic.phasePoint0() == 1 && mic.coordPoint0(axis) != mic.axisSize(axis) - 1) {
//                        lmax++;
//                        mic.advancePoint0(axis);
//                    }
//
//                    for (int r = 0; r <= lmax; r++) nonPeriodic[r]++;
//
//                }
//            }
//        }
//
//        for (int r = 0; r < rmax; r++) {
//            nonPeriodic[r] = nonPeriodic[r] / den[r];
//        }
//
//        return nonPeriodic;
//    }
//
//};
//
//class TwoPoint : public CorrFunction {
//
//public:
//    vector<int> aux;
//    vector<int> curr;
//    vector<double> target;
//
//    Microstructure& mic;
//    int rmax, numAxis;
//    double currEnergy, auxEnergy;
//    bool isReference = false;
//
//
//    TwoPoint(Microstructure& mic_, int numAxis_) : mic(mic_) { //assume como referência
//        isReference = true;
//
//        numAxis = numAxis_;
//
//        if (numAxis == 1) {
//            rmax = mic.nx / 2;
//        }
//        else if (numAxis == 2) {
//            rmax = min(mic.nx / 2, mic.ny / 2);
//        }
//        else if (numAxis == 3) {
//            rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
//        }
//
//        firstCalc();
//
//        target.assign(rmax, 0);
//
//        for (int r = 0; r < rmax; r++) {
//            target[r] = double(curr[r]) / (numAxis * mic.n);
//        }
//
//
//        currEnergy = getStartEnergy();
//
//    }
//
//    TwoPoint(Microstructure& mic_, vector<double> target_, int numAxis_) : mic(mic_) { //lê do vetor
//        numAxis = numAxis_;
//
//        if (numAxis == 1) {
//            rmax = mic.nx / 2;
//        }
//        else if (numAxis == 2) {
//            rmax = min(mic.nx / 2, mic.ny / 2);
//        }
//        else if (numAxis == 3) {
//            rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
//        }
//
//
//        firstCalc();
//        target = target_;
//        currEnergy = getStartEnergy();
//
//    }
//
//
//
//
//
//    TwoPoint(Microstructure& mic_, string targetPath, int numAxis_) : mic(mic_) {
//        numAxis = numAxis_;
//
//        if (numAxis == 1) {
//            rmax = mic.nx / 2;
//        }
//        else if (numAxis == 2) {
//            rmax = min(mic.nx / 2, mic.ny / 2);
//        }
//        else if (numAxis == 3) {
//            rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
//        }
//
//
//        firstCalc();
//        readFile(targetPath);
//        currEnergy = getStartEnergy();
//
//    }
//
//    void firstCalc() {
//
//        curr.assign(rmax, 0);
//
//        for (int axis = 0; axis < numAxis; axis++) {
//            for (int i = 0; i < mic.n; i++) {
//                mic.setPoint0(i);
//                mic.setPoint1(i);
//                for (int r = 0; r < rmax; r++) {
//                    if (mic.phasePoint0() == 1 and mic.phasePoint1() == 1) {
//                        curr[r]++;
//                    }
//                    mic.advancePoint1(axis);
//                }
//            }
//        }
//    }
//
//    double getStartEnergy() {
//        double sum = 0;
//        for (int r = 0; r < rmax; r++) {
//            sum += pow(target[r] - double(curr[r]) / (numAxis * mic.n), 2);
//        }
//        return sum;
//    }
//
//    double getAuxEnergy() override {
//        double sum = 0;
//        for (int r = 0; r < rmax; r++) {
//            sum += pow(target[r] - double(aux[r]) / (numAxis * mic.n), 2);
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
//        string value;
//
//        string line;
//
//        target.assign(rmax, 0);
//
//        for (int r = 0; r < rmax; r++) {
//            file >> value >> value;
//            target[r] = stod(value);
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
//        //troca
//        mic.arr[n0_] = 1;
//        mic.arr[n1_] = 0;
//
//
//        for (int axis = 0; axis < numAxis; axis++) {
//
//            mic.setPoint0(n0_);
//            mic.setPoint1(n0_);
//
//            for (int r = 1; r < rmax; r++) {
//                mic.backPoint0(axis);
//                mic.advancePoint1(axis);
//
//                if (mic.phasePoint0() == 1) aux[r]++;
//                if (mic.phasePoint1() == 1) aux[r]++;
//            }
//
//            mic.setPoint0(n1_);
//            mic.setPoint1(n1_);
//
//            for (int r = 1; r < rmax; r++) {
//                mic.backPoint0(axis);
//                mic.advancePoint1(axis);
//
//                if (mic.phasePoint0() == 1) aux[r]--;
//                if (mic.phasePoint1() == 1) aux[r]--;
//            }
//
//            //destroca       
//            mic.arr[n0_] = 0;
//            mic.arr[n1_] = 1;
//
//            auxEnergy = getAuxEnergy();
//
//
//        }
//    }
//
//    double getCurrEnergy() override {
//        return currEnergy;
//    }
//
//    void writeToFile(string path_) {
//        ofstream file;
//        double den = numAxis * mic.n;
//        file.open(path_);
//        int pos = 0;
//
//        for (int r = 0; r < rmax; r++)
//        {
//            file << r << "\t" << curr[r] / den << "\n";
//            pos++;
//        }
//        file.close();
//    }
//
//    void writeToFileNonPeriodic(string path_) {
//        ofstream file;
//        vector<double> curr = evalNonPeriodic();
//        file.open(path_);
//        int pos = 0;
//
//        for (int r = 0; r < rmax; r++)
//        {
//            file << r << "\t" << curr[r] << "\n";
//            pos++;
//        }
//        file.close();
//    }
//
//    vector<double> evalNonPeriodic() {
//
//        vector<double> nonPeriodic;
//        vector<double> den;
//        nonPeriodic.assign(rmax, 0);
//        den.assign(rmax, 0);
//
//        int diff = 0;
//
//        for (int axis = 0; axis < numAxis; axis++) {
//
//            for (int i = 0; i < mic.n; i++) {
//                mic.setPoint0(i);
//                mic.setPoint1(i);
//                diff = 0;
//                while (diff < rmax) {
//                    if (mic.phasePoint0() == 1 and mic.phasePoint1() == 1) nonPeriodic[diff]++;
//
//                    mic.advancePoint1(axis);
//                    den[diff]++;
//                    diff++;
//
//                    if (mic.coordPoint1(axis) == mic.axisSize(axis) - 1) break;
//
//                }
//            }
//        }
//
//        for (int r = 0; r < rmax; r++) {
//            nonPeriodic[r] = nonPeriodic[r] / den[r];
//        }
//
//        return nonPeriodic;
//    }
//
//};
