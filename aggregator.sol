// Solidity program to
// demonstrate retrieve
// values from the mapping
pragma solidity >=0.6.12;



contract ExperimentAgreegator {

    uint256 constant private U255_MAX_PLUS_1 = 57896044618658097711785492504343953926634992332820282019728792003956564819968;

    /// @dev Modular euclidean inverse of a number (mod p).
    /// @param _x The number
    /// @param _pp The modulus
    /// @return q such that x*q = 1 (mod _pp)
    function invMod(uint256 _x, uint256 _pp) internal pure returns (uint256) {
      require(_x != 0 && _x != _pp && _pp != 0, "Invalid number");
      uint256 q = 0;
      uint256 newT = 1;
      uint256 r = _pp;
      uint256 t;
      while (_x != 0) {
        t = r / _x;
        (q, newT) = (newT, addmod(q, (_pp - mulmod(t, newT, _pp)), _pp));
        (r, _x) = (_x, r - t * _x);
      }
  
      return q;
    }
  
    /// @dev Modular exponentiation, b^e % _pp.
    /// Source: https://github.com/androlo/standard-contracts/blob/master/contracts/src/crypto/ECCMath.sol
    /// @param _base base
    /// @param _exp exponent
    /// @param _pp modulus
    /// @return r such that r = b**e (mod _pp)
    function expMod(uint256 _base, uint256 _exp, uint256 _pp) internal pure returns (uint256) {
      require(_pp!=0, "Modulus is zero");
  
      if (_base == 0)
        return 0;
      if (_exp == 0)
        return 1;
  
      uint256 r = 1;
      uint256 bit = U255_MAX_PLUS_1;
      assembly {
        for { } gt(bit, 0) { }{
          r := mulmod(mulmod(r, r, _pp), exp(_base, iszero(iszero(and(_exp, bit)))), _pp)
          r := mulmod(mulmod(r, r, _pp), exp(_base, iszero(iszero(and(_exp, div(bit, 2))))), _pp)
          r := mulmod(mulmod(r, r, _pp), exp(_base, iszero(iszero(and(_exp, div(bit, 4))))), _pp)
          r := mulmod(mulmod(r, r, _pp), exp(_base, iszero(iszero(and(_exp, div(bit, 8))))), _pp)
          bit := div(bit, 16)
        }
      }
  
      return r;
    }
  
    /// @dev Converts a point (x, y, z) expressed in Jacobian coordinates to affine coordinates (x', y', 1).
    /// @param _x coordinate x
    /// @param _y coordinate y
    /// @param _z coordinate z
    /// @param _pp the modulus
    /// @return (x', y') affine coordinates
    function toAffine(
      uint256 _x,
      uint256 _y,
      uint256 _z,
      uint256 _pp)
    internal pure returns (uint256, uint256)
    {
      uint256 zInv = invMod(_z, _pp);
      uint256 zInv2 = mulmod(zInv, zInv, _pp);
      uint256 x2 = mulmod(_x, zInv2, _pp);
      uint256 y2 = mulmod(_y, mulmod(zInv, zInv2, _pp), _pp);
  
      return (x2, y2);
    }
  
    /// @dev Derives the y coordinate from a compressed-format point x [[SEC-1]](https://www.secg.org/SEC1-Ver-1.0.pdf).
    /// @param _prefix parity byte (0x02 even, 0x03 odd)
    /// @param _x coordinate x
    /// @param _aa constant of curve
    /// @param _bb constant of curve
    /// @param _pp the modulus
    /// @return y coordinate y
    function deriveY(
      uint8 _prefix,
      uint256 _x,
      uint256 _aa,
      uint256 _bb,
      uint256 _pp)
    internal pure returns (uint256)
    {
      require(_prefix == 0x02 || _prefix == 0x03, "Invalid compressed EC point prefix");
  
      // x^3 + ax + b
      uint256 y2 = addmod(mulmod(_x, mulmod(_x, _x, _pp), _pp), addmod(mulmod(_x, _aa, _pp), _bb, _pp), _pp);
      y2 = expMod(y2, (_pp + 1) / 4, _pp);
      // uint256 cmp = yBit ^ y_ & 1;
      uint256 y = (y2 + _prefix) % 2 == 0 ? y2 : _pp - y2;
  
      return y;
    }
  
    /// @dev Check whether point (x,y) is on curve defined by a, b, and _pp.
    /// @param _x coordinate x of P1
    /// @param _y coordinate y of P1
    /// @param _aa constant of curve
    /// @param _bb constant of curve
    /// @param _pp the modulus
    /// @return true if x,y in the curve, false else
    function isOnCurve(
      uint _x,
      uint _y,
      uint _aa,
      uint _bb,
      uint _pp)
    internal pure returns (bool)
    {
      if (0 == _x || _x >= _pp || 0 == _y || _y >= _pp) {
        return false;
      }
      // y^2
      uint lhs = mulmod(_y, _y, _pp);
      // x^3
      uint rhs = mulmod(mulmod(_x, _x, _pp), _x, _pp);
      if (_aa != 0) {
        // x^3 + a*x
        rhs = addmod(rhs, mulmod(_x, _aa, _pp), _pp);
      }
      if (_bb != 0) {
        // x^3 + a*x + b
        rhs = addmod(rhs, _bb, _pp);
      }
  
      return lhs == rhs;
    }
  
    /// @dev Calculate inverse (x, -y) of point (x, y).
    /// @param _x coordinate x of P1
    /// @param _y coordinate y of P1
    /// @param _pp the modulus
    /// @return (x, -y)
    function ecInv(
      uint256 _x,
      uint256 _y,
      uint256 _pp)
    internal pure returns (uint256, uint256)
    {
      return (_x, (_pp - _y) % _pp);
    }
  
    /// @dev Add two points (x1, y1) and (x2, y2) in affine coordinates.
    /// @param _x1 coordinate x of P1
    /// @param _y1 coordinate y of P1
    /// @param _x2 coordinate x of P2
    /// @param _y2 coordinate y of P2
    /// @param _aa constant of the curve
    /// @param _pp the modulus
    /// @return (qx, qy) = P1+P2 in affine coordinates
    function ecAdd(
      uint256 _x1,
      uint256 _y1,
      uint256 _x2,
      uint256 _y2,
      uint256 _aa,
      uint256 _pp)
      internal pure returns(uint256, uint256)
    {
      uint x = 0;
      uint y = 0;
      uint z = 0;
  
      // Double if x1==x2 else add
      if (_x1==_x2) {
        // y1 = -y2 mod p
        if (addmod(_y1, _y2, _pp) == 0) {
          return(0, 0);
        } else {
          // P1 = P2
          (x, y, z) = jacDouble(
            _x1,
            _y1,
            1,
            _aa,
            _pp);
        }
      } else {
        (x, y, z) = jacAdd(
          _x1,
          _y1,
          1,
          _x2,
          _y2,
          1,
          _pp);
      }
      // Get back to affine
      return toAffine(
        x,
        y,
        z,
        _pp);
    }
  
    /// @dev Substract two points (x1, y1) and (x2, y2) in affine coordinates.
    /// @param _x1 coordinate x of P1
    /// @param _y1 coordinate y of P1
    /// @param _x2 coordinate x of P2
    /// @param _y2 coordinate y of P2
    /// @param _aa constant of the curve
    /// @param _pp the modulus
    /// @return (qx, qy) = P1-P2 in affine coordinates
    function ecSub(
      uint256 _x1,
      uint256 _y1,
      uint256 _x2,
      uint256 _y2,
      uint256 _aa,
      uint256 _pp)
    internal pure returns(uint256, uint256)
    {
      // invert square
      (uint256 x, uint256 y) = ecInv(_x2, _y2, _pp);
      // P1-square
      return ecAdd(
        _x1,
        _y1,
        x,
        y,
        _aa,
        _pp);
    }
  
    /// @dev Multiply point (x1, y1, z1) times d in affine coordinates.
    /// @param _k scalar to multiply
    /// @param _x coordinate x of P1
    /// @param _y coordinate y of P1
    /// @param _aa constant of the curve
    /// @param _pp the modulus
    /// @return (qx, qy) = d*P in affine coordinates
    function ecMul(
      uint256 _k,
      uint256 _x,
      uint256 _y,
      uint256 _aa,
      uint256 _pp)
    internal pure returns(uint256, uint256)
    {
      // Jacobian multiplication
      (uint256 x1, uint256 y1, uint256 z1) = jacMul(
        _k,
        _x,
        _y,
        1,
        _aa,
        _pp);
      // Get back to affine
      return toAffine(
        x1,
        y1,
        z1,
        _pp);
    }
  
    /// @dev Adds two points (x1, y1, z1) and (x2 y2, z2).
    /// @param _x1 coordinate x of P1
    /// @param _y1 coordinate y of P1
    /// @param _z1 coordinate z of P1
    /// @param _x2 coordinate x of square
    /// @param _y2 coordinate y of square
    /// @param _z2 coordinate z of square
    /// @param _pp the modulus
    /// @return (qx, qy, qz) P1+square in Jacobian
    function jacAdd(
      uint256 _x1,
      uint256 _y1,
      uint256 _z1,
      uint256 _x2,
      uint256 _y2,
      uint256 _z2,
      uint256 _pp)
    internal pure returns (uint256, uint256, uint256)
    {
      if (_x1==0 && _y1==0)
        return (_x2, _y2, _z2);
      if (_x2==0 && _y2==0)
        return (_x1, _y1, _z1);
  
      // We follow the equations described in https://pdfs.semanticscholar.org/5c64/29952e08025a9649c2b0ba32518e9a7fb5c2.pdf Section 5
      uint[4] memory zs; // z1^2, z1^3, z2^2, z2^3
      zs[0] = mulmod(_z1, _z1, _pp);
      zs[1] = mulmod(_z1, zs[0], _pp);
      zs[2] = mulmod(_z2, _z2, _pp);
      zs[3] = mulmod(_z2, zs[2], _pp);
  
      // u1, s1, u2, s2
      zs = [
        mulmod(_x1, zs[2], _pp),
        mulmod(_y1, zs[3], _pp),
        mulmod(_x2, zs[0], _pp),
        mulmod(_y2, zs[1], _pp)
      ];
  
      // In case of zs[0] == zs[2] && zs[1] == zs[3], double function should be used
      require(zs[0] != zs[2] || zs[1] != zs[3], "Use jacDouble function instead");
  
      uint[4] memory hr;
      //h
      hr[0] = addmod(zs[2], _pp - zs[0], _pp);
      //r
      hr[1] = addmod(zs[3], _pp - zs[1], _pp);
      //h^2
      hr[2] = mulmod(hr[0], hr[0], _pp);
      // h^3
      hr[3] = mulmod(hr[2], hr[0], _pp);
      // qx = -h^3  -2u1h^2+r^2
      uint256 qx = addmod(mulmod(hr[1], hr[1], _pp), _pp - hr[3], _pp);
      qx = addmod(qx, _pp - mulmod(2, mulmod(zs[0], hr[2], _pp), _pp), _pp);
      // qy = -s1*z1*h^3+r(u1*h^2 -x^3)
      uint256 qy = mulmod(hr[1], addmod(mulmod(zs[0], hr[2], _pp), _pp - qx, _pp), _pp);
      qy = addmod(qy, _pp - mulmod(zs[1], hr[3], _pp), _pp);
      // qz = h*z1*z2
      uint256 qz = mulmod(hr[0], mulmod(_z1, _z2, _pp), _pp);
      return(qx, qy, qz);
    }
  
    /// @dev Doubles a points (x, y, z).
    /// @param _x coordinate x of P1
    /// @param _y coordinate y of P1
    /// @param _z coordinate z of P1
    /// @param _aa the a scalar in the curve equation
    /// @param _pp the modulus
    /// @return (qx, qy, qz) 2P in Jacobian
    function jacDouble(
      uint256 _x,
      uint256 _y,
      uint256 _z,
      uint256 _aa,
      uint256 _pp)
    internal pure returns (uint256, uint256, uint256)
    {
      if (_z == 0)
        return (_x, _y, _z);
  
      // We follow the equations described in https://pdfs.semanticscholar.org/5c64/29952e08025a9649c2b0ba32518e9a7fb5c2.pdf Section 5
      // Note: there is a bug in the paper regarding the m parameter, M=3*(x1^2)+a*(z1^4)
      // x, y, z at this point represent the squares of _x, _y, _z
      uint256 x = mulmod(_x, _x, _pp); //x1^2
      uint256 y = mulmod(_y, _y, _pp); //y1^2
      uint256 z = mulmod(_z, _z, _pp); //z1^2
  
      // s
      uint s = mulmod(4, mulmod(_x, y, _pp), _pp);
      // m
      uint m = addmod(mulmod(3, x, _pp), mulmod(_aa, mulmod(z, z, _pp), _pp), _pp);
  
      // x, y, z at this point will be reassigned and rather represent qx, qy, qz from the paper
      // This allows to reduce the gas cost and stack footprint of the algorithm
      // qx
      x = addmod(mulmod(m, m, _pp), _pp - addmod(s, s, _pp), _pp);
      // qy = -8*y1^4 + M(S-T)
      y = addmod(mulmod(m, addmod(s, _pp - x, _pp), _pp), _pp - mulmod(8, mulmod(y, y, _pp), _pp), _pp);
      // qz = 2*y1*z1
      z = mulmod(2, mulmod(_y, _z, _pp), _pp);
  
      return (x, y, z);
    }
  
    /// @dev Multiply point (x, y, z) times d.
    /// @param _d scalar to multiply
    /// @param _x coordinate x of P1
    /// @param _y coordinate y of P1
    /// @param _z coordinate z of P1
    /// @param _aa constant of curve
    /// @param _pp the modulus
    /// @return (qx, qy, qz) d*P1 in Jacobian
    function jacMul(
      uint256 _d,
      uint256 _x,
      uint256 _y,
      uint256 _z,
      uint256 _aa,
      uint256 _pp)
    internal pure returns (uint256, uint256, uint256)
    {
      // Early return in case that `_d == 0`
      if (_d == 0) {
        return (_x, _y, _z);
      }
  
      uint256 remaining = _d;
      uint256 qx = 0;
      uint256 qy = 0;
      uint256 qz = 1;
  
      // Double and add algorithm
      while (remaining != 0) {
        if ((remaining & 1) != 0) {
          (qx, qy, qz) = jacAdd(
            qx,
            qy,
            qz,
            _x,
            _y,
            _z,
            _pp);
        }
        remaining = remaining / 2;
        (_x, _y, _z) = jacDouble(
          _x,
          _y,
          _z,
          _aa,
          _pp);
      }
      return (qx, qy, qz);
    }


///---------------------Contract's Memory---------------------\\\

//Weight related stuff
struct weight{
    uint times;
    uint id;
}

mapping (uint => weight) barh;

function initweights(uint n, uint timi) public payable{
    for(uint i=0; i<n; i++) { barh[i].times=timi;}
}

uint256 public howmany=0;

// function setWeightsAdditively(uint256 posa, uint256[] tabarh) public payable {
//     for(uint i=howmany; i<howmany+posa; i++){ barh[i].times = tabarh[i]; }howmany = howmany + posa;
// }

function changeSpecigicWeight(uint256 index, uint256 value) public payable {
    barh[index].times = value;
}

function getWeights(uint256 index/*, uint256 adID*/) public view returns(uint256 baros){
    //require(devices.address[adID]==msg.sender);
    return barh[index].times;
}

function resetWeights(uint n) public payable{
     for(uint i=0; i<n; i++) { delete(barh[i].times);}
}

// Master secret key
//uint256 msk = 0;
uint256 public constant aa = 0;
// Field size
uint256 constant pp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;

// Base point (generator) G
uint256[2] G = [0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8];

// Order of G
uint256 constant nn = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;

// Maximum value of s
uint256 constant lowSmax = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5D576E7357A4501DDFE92F46681B20A0;


// List of timers that each phase MUST end by an explicit time in UNIX timestamp.
// Ethereum works in SECONDS. Not milliseconds.
//uint  finishSETUPPhase; // SC to transition to next phase.
uint256 endSETUPPhase; // SC does not transition to next phase by this time.
uint256 endCIPHERCOLLECTIONPhase; // Devices have not sent their ciphertexts in by this time.
uint256 endCIPHERCSUMATIONPhase; // Devices have not submitted their ciphertexts by this stage.
uint256 gap = 500;

uint256 public totaleligible;
uint256 public totalsubmitted;
uint256 public totalremaining;
uint256 totalrequested;

uint256[2] PKD;

///---------------------Struct---------------------\\\

struct Device   {
    uint id;
    address andres;
    uint256[3] PUBKEY;
    uint256[2] CC;
    uint256[2] CD;
    bool eligible;
    bool ciphercast;
    }

mapping (uint => Device) public devices;
mapping (address => uint) public addressid; // Address to Counter

///---------------------States & Modifiers---------------------\\\

enum  State { SETUP, CIPHERCOLLECTION, CIPHERCSUMATION, FINISHED }
State public state;

modifier inState(State s) {
    if(state != s) {
        revert();
        }
    _;
    }

modifier eligibleforever() {
    bool eli = false;
    for( uint i=0; i<totaleligible; i++){
        if (msg.sender == devices[i].andres) { eli = true; }
        }
    if (eli == false) { revert(); }
    _;
    }

/*    address addrDecryptor;
modifier onlyDecryptor {
    if(addrDecryptor != msg.sender) revert();
    _;
    }
*/
///---------------------Functions---------------------\\\

function setEligible(address addr/*, uint256[3] _PK*/)  inState(State.SETUP) public payable returns (uint256 a) {
    // Sign up the addresses
    //for(uint i=totaleligible; i<totaleligible+addr.length; i++) {
        if(!devices[totaleligible].eligible) {
            addressid[addr] = totaleligible;
            devices[totaleligible].andres = addr;
            //devices[totaleligible].PUBKEY = _PK;
            devices[totaleligible].eligible = true;
            totaleligible += 1;
            }
     //   }
    return totaleligible;
    }

// Owner of contract declares that eligible addresses begin round 1 of the protocol
// Time is the number of 'blocks' we must wait until we can move onto round 2.
uint tx_cost;
function beginProc(uint _endCIPHERCOLLECTIONPhase, uint _endCIPHERCSUMATIONPhase, uint budget, uint costPerTransaction) inState(State.SETUP)   public payable returns (bool){
    require(msg.value == budget);
    state = State.CIPHERCOLLECTION;
    endCIPHERCOLLECTIONPhase = block.timestamp + _endCIPHERCOLLECTIONPhase;
    endCIPHERCSUMATIONPhase = block.timestamp + _endCIPHERCSUMATIONPhase;
    tx_cost = costPerTransaction;
    return true;
    }

function collectionDeadlinePassed() public returns (bool){
    if(state == State.CIPHERCOLLECTION && block.timestamp > endCIPHERCOLLECTIONPhase) {
        // Check which devices have not given their ciphertexts
        for(uint i=0; i<totaleligible; i++) {
            if(devices[i].ciphercast) {
                totalremaining++;
                }
            }
        if(totalremaining == 0){
            state = State.CIPHERCSUMATION;
            prodCipher();
            }

         /*   else{
        //Wait for randomly elected device to send encrypted zeros
        //Maybe have devices send as well encrypted zeros and store them on chain and use these values to fill the gaps IS IT REALLY NECESSARY?
        state = State.FINISHED;
        }
        return true; */ //For now we assume that everyone sends their encrypted input
        }
    }

function resetToSetup() public   {
    state = State.SETUP;
    for (uint h=0; h<totaleligible; h++){
        delete(devices[h].CC);
        devices[h].ciphercast = false;
        devices[h].eligible = false;
        delete(totalsubmitted);
        }
    delete(totaleligible);
    }

function goToCipherCollection() public   {
    state = State.CIPHERCOLLECTION;
    }

function goToCipherSumation() public   {
    state = State.CIPHERCSUMATION;
//    prodCipher();
    }

//---------------------CRYPTO---------------------//



// Change andres to msg.sender
uint256[3] random;
uint256 public roundcounter=0;

uint public counterdg =0;
function resetcounterdg() public {
counterdg =0;
}
function simpleappendCiphers(uint256[2] memory appendMe/*, uint256[2] appendalsoMe*/)  inState(State.CIPHERCOLLECTION)/*eligibleforever()*/ public payable {
    devices[counterdg].ciphercast = true;
    devices[counterdg].CC=appendMe;
    totalsubmitted += 1;
    counterdg +=1;
    }
function appendCiphers(address ad, uint256[2] memory appendMe/*, uint256[2] appendalsoMe*/)  inState(State.CIPHERCOLLECTION)/*eligibleforever()*/ public payable {
    devices[counterdg].ciphercast = true;
    devices[counterdg].CC=appendMe;
    totalsubmitted += 1;
    counterdg +=1;
    if (totalsubmitted == totaleligible){
        goToCipherSumation();
        }
    payable(msg.sender).transfer(tx_cost);
    }


uint256[2] testpoint;
function testAddG() inState(State.CIPHERCSUMATION) public payable returns ( uint256, uint256){
    ( testpoint[0] ,testpoint[1])=ecAdd(G[0], G[1], G[0], G[1], aa, pp);
    return (testpoint[0] ,testpoint[1]);
}
function testmulG() inState(State.CIPHERCSUMATION) public payable  returns ( uint256, uint256){
    (testpoint[0] ,testpoint[1])=ecMul(10000, G[0],  G[1], aa, pp);
    return (testpoint[0] ,testpoint[1]);
}
 function gettestpoint() public view returns(uint256, uint256){
    //require(devices.address[adID]==msg.sender);
    return (testpoint[0] ,testpoint[1]);
    }



uint256[3] public FINALCIPHERTEXTC;
function prodCipher() inState(State.CIPHERCSUMATION) public payable returns ( uint256, uint256) {  
    (FINALCIPHERTEXTC[0],FINALCIPHERTEXTC[1])=ecMul(barh[0].times, devices[0].CC[0], devices[0].CC[1], aa, pp);
    for (uint i=1; i<totaleligible;  i++) {
        ( uint256 x, uint256 y)=ecMul(barh[i].times, devices[i].CC[0], devices[i].CC[1], aa, pp);
        (FINALCIPHERTEXTC[0],FINALCIPHERTEXTC[1]) = ecAdd(FINALCIPHERTEXTC[0], FINALCIPHERTEXTC[1], x, y, aa, pp);
        }
    return (FINALCIPHERTEXTC[0],FINALCIPHERTEXTC[1]);
    }

//totaleligible
// function prodCipherWithParams(uint apo, uint ews) inState(State.CIPHERCSUMATION) public payable returns (uint256[3] FC) {
//     FINALCIPHERTEXTC = devices[0].CC;
//     for (uint i=apo; i<ews; i++) {
//         FINALCIPHERTEXTC = _add(FINALCIPHERTEXTC,_mul(barh[i].times,devices[i].CC));
//         }
//     toZ1(FINALCIPHERTEXTC,pp);
//     return (FINALCIPHERTEXTC);
//     }



// function getDevice(uint256 index)   public view returns (uint[3] _registeredkey, uint[2] _ciphertextc, uint[2] _ciphertextd, address _address){
//     // uint index = addressid[adreas];
//     _registeredkey = devices[index].PUBKEY;
//     _ciphertextc = devices[index].CC;
//     _ciphertextd = devices[index].CD;
//     _address = devices[index].andres;
//     }


function getfCiphertext() public view returns(uint[3] memory FCIPHERC){
    return (FINALCIPHERTEXTC);
    }


// function getgtotheri() public view returns(uint256[3] GR){
//     return random;
//     }


// function setHowManyEligible(uint256 hm) public payable{
//     totaleligible=hm;
// }


// uint256[2] gS;
// uint256[2] MPK;
// uint256[3] MTOTAL;
// uint256[2] gA1S;


// function isEqual(uint256[3] A, uint256[3] B) public payable returns (bool t){
//     if ((A[0] == B[0]) && (A[1] == B[1]) && (A[1] == B[1])) { return true; }
// }

}