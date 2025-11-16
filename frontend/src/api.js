// src/api.js
export async function runBeamSimulation(params) {
    const response = await fetch("http://localhost:8000/api/run_beam", {  // or full http://localhost:8000/api/run_beam if no proxy
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(params),
    });

    let errorDetail = "";
    if (!response.ok) {
        try {
            const errJson = await response.json();
            if (errJson.detail) errorDetail = `: ${errJson.detail}`;
        } catch { }
        throw new Error(`Backend error: ${response.status}${errorDetail}`);
    }
    return await response.json();
}