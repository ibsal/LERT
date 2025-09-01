// ---- Simulated hydrostatic pressure model with serial control ----

float pressure_psi = 0.0f;        // current pressure (psi)
int   valve_open_pct = 0;         // 0..100
const float P_MAX = 1500.0f;      // clamp (for safety)

const float FILL_PSI_PER_SEC = 4.0f;   // closed fill rate (psi/s) at 0% open
const float K_MAX = 6.0f;              // max vent coefficient (1/s) at 100% open
const float NOISE_PSI = 0.02f;         // small noise amplitude (psi)

unsigned long last_update_ms = 0;

void setup() {
  Serial.begin(9600);
  while (!Serial) { /* wait for port */ }
  last_update_ms = millis();
}

void updatePressure() {
  unsigned long now = millis();
  float dt = (now - last_update_ms) * 1e-3f;   // seconds
  if (dt <= 0) return;
  last_update_ms = now;

  // normalized openness [0..1]
  float u = constrain(valve_open_pct / 100.0f, 0.0f, 1.0f);

  // Fill term: slower as valve opens (venting opposes filling)
  float dP_fill = FILL_PSI_PER_SEC * (1.0f - u);

  // Vent term: proportional to openness and current pressure (exponential decay)
  float dP_vent = - K_MAX * u * pressure_psi;

  // Integrate
  pressure_psi += (dP_fill + dP_vent) * dt;

  // Add a touch of noise
  pressure_psi += NOISE_PSI * (float)random(-5, 6) / 5.0f;

  // Clamp
  if (pressure_psi < 0.0f) pressure_psi = 0.0f;
  if (pressure_psi > P_MAX) pressure_psi = P_MAX;
}

// Helper to read a trimmed line (blocking briefly until newline)
String readLineTrimmed() {
  String s = Serial.readStringUntil('\n');
  s.trim();
  return s;
}

void loop() {
  // Always advance the simulation
  updatePressure();

  // Handle commands if any
  if (Serial.available()) {
    String cmd = readLineTrimmed();

    if (cmd == "COMCHECK") {
      Serial.println("HERE");

    } else if (cmd == "PRESSURE") {
      Serial.println(pressure_psi, 4);

    } else if (cmd == "TIME") {
      Serial.println(millis());

    } else if (cmd == "CLOSE") {
      // Fully closed: slow build-up
      valve_open_pct = 0;

    } else if (cmd == "OPEN") {
      // Next line should be an integer 0..100 sent by MATLAB
      // Wait briefly if needed for the next line to arrive
      unsigned long t0 = millis();
      while (!Serial.available() && (millis() - t0) < 100) { /* short wait */ }
      String pctLine = readLineTrimmed();
      int pct = pctLine.toInt();          // toInt() handles non-numeric gracefully (returns 0)
      valve_open_pct = constrain(pct, 0, 100);

    } else {
      // Unknown command: ignore
    }
  }
}
