#include <Servo.h>

// ========= USER SETTINGS =========
// Pressure transducer (same as before)
const uint8_t  SENSOR_PIN = A0;
const float V_MIN = 0.5f;
const float V_MAX = 4.5f;
const float PSI_FULL_SCALE = 1000.0f;
const int   ADC_SAMPLES = 16;
const float EMA_ALPHA   = 0.25f;

// Servo pin
const uint8_t SERVO_PIN = 2;

// ---- Mechanical servo endpoints in degrees (can be swapped) ----
// These are the ACTUAL angles that hit your valve stops.
// If 26° is fully CLOSED and 104° is fully OPEN, set like below.
// If opposite, just swap them (works either way).
int SERVO_CLOSED_DEG = 26;
int SERVO_OPEN_DEG   = 104;

// Optional: define the MATLAB command range that should map to closed/open.
// If MATLAB always sends 0..100, leave these as 0 and 100.
int CMD_PCT_MIN = 0;    // percent that maps to SERVO_CLOSED_DEG
int CMD_PCT_MAX = 100;  // percent that maps to SERVO_OPEN_DEG
// =================================

Servo valve;
int   last_cmd_pct = 0;
float ema_pressure_psi = 0.0f;
unsigned long last_update_ms = 0;

// -------- PT helpers --------
static inline float readVoltage() {
  long acc = 0;
  for (int i = 0; i < ADC_SAMPLES; ++i) acc += analogRead(SENSOR_PIN);
  float counts = (float)acc / ADC_SAMPLES;
  return counts * (5.0f / 1023.0f);
}

float readPressurePsi() {
  float v = readVoltage();
  float span = (V_MAX - V_MIN);
  float psi = 0.0f;

  if (span > 0.001f) {
    float x = (v - V_MIN) / span;   // expected 0..1
    if (x < 0.0f) x = 0.0f;
    if (x > 1.0f) x = 1.0f;
    psi = x * PSI_FULL_SCALE;
  }

  static bool first = true;
  if (first) { ema_pressure_psi = psi; first = false; }
  else       { ema_pressure_psi = EMA_ALPHA * psi + (1.0f - EMA_ALPHA) * ema_pressure_psi; }
  return ema_pressure_psi;
}

// -------- Valve mapping --------
// Normalize MATLAB percent to 0..1 using CMD_PCT_MIN/MAX (handles swapped)
float pctTo01(int pct) {
  if (CMD_PCT_MAX == CMD_PCT_MIN) return 0.0f;
  float f = (float)(pct - CMD_PCT_MIN) / (float)(CMD_PCT_MAX - CMD_PCT_MIN);
  if (f < 0.0f) f = 0.0f;
  if (f > 1.0f) f = 1.0f;
  return f;
}

// Map 0..1 to servo angle between the two endpoints (handles swapped)
int map01ToAngle(float f01) {
  int a0 = constrain(SERVO_CLOSED_DEG, 0, 180);
  int a1 = constrain(SERVO_OPEN_DEG,   0, 180);
  float angle = a0 + f01 * (float)(a1 - a0);
  int adeg = (int)(angle + (angle >= 0 ? 0.5f : -0.5f));
  return constrain(adeg, 0, 180);
}

void applyValvePercent(int pct) {
  last_cmd_pct = pct;
  float f01 = pctTo01(pct);
  int angle = map01ToAngle(f01);
  valve.write(angle);  // degrees
}

// -------- Serial helpers --------
String readLineTrimmed() {
  String s = Serial.readStringUntil('\n');
  s.trim();
  return s;
}

// -------- Setup / Loop --------
void setup() {
  Serial.begin(9600);
  while (!Serial) { /* wait for port */ }

  analogReference(DEFAULT);
  pinMode(SENSOR_PIN, INPUT);

  valve.attach(SERVO_PIN);
  // Go to fully "closed" mechanical angle at startup:
  valve.write(constrain(SERVO_OPEN_DEG, 0, 180));

  last_update_ms = millis();
}

void loop() {
  if (Serial.available()) {
    String cmd = readLineTrimmed();

    if (cmd == "COMCHECK") {
      Serial.println("HERE");

    } else if (cmd == "PRESSURE") {
      Serial.println(readPressurePsi(), 4);

    } else if (cmd == "TIME") {
      Serial.println(millis());

    } else if (cmd == "CLOSE") {
      // Force fully-closed mechanical angle
      valve.write(constrain(SERVO_CLOSED_DEG, 0, 180));

    } else if (cmd == "OPEN") {
      // Next line = MATLAB "percent" (any integer); map to 26°..104° range
      unsigned long t0 = millis();
      while (!Serial.available() && (millis() - t0) < 150) { /* wait a bit */ }
      int pct = readLineTrimmed().toInt();
      applyValvePercent(pct);

    // --- Optional runtime calibration commands ---
    } else if (cmd == "SET_CLOSE_ANGLE") {
      SERVO_CLOSED_DEG = constrain(readLineTrimmed().toInt(), 0, 180);
      applyValvePercent(last_cmd_pct);

    } else if (cmd == "SET_OPEN_ANGLE") {
      SERVO_OPEN_DEG = constrain(readLineTrimmed().toInt(), 0, 180);
      applyValvePercent(last_cmd_pct);

    } else if (cmd == "SET_CMD_MIN") {  // percent that maps to closed angle
      CMD_PCT_MIN = readLineTrimmed().toInt();
      applyValvePercent(last_cmd_pct);

    } else if (cmd == "SET_CMD_MAX") {  // percent that maps to open angle
      CMD_PCT_MAX = readLineTrimmed().toInt();
      applyValvePercent(last_cmd_pct);
    }
  }

  // keep pressure filter warm
  if (millis() - last_update_ms >= 50) {
    readPressurePsi();
    last_update_ms = millis();
  }
}
