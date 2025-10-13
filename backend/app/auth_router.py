from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel
import os
from google.oauth2 import id_token
from google.auth.transport import requests

router = APIRouter()


class GoogleLoginRequest(BaseModel):
    token: str


# Optional audience verification for additional security. Set GOOGLE_CLIENT_ID in env.
GOOGLE_CLIENT_ID = os.getenv("GOOGLE_CLIENT_ID")

@router.post("/google-login")
async def google_login(payload: GoogleLoginRequest):
    """
    Verify Google login token and return basic user info
    """
    try:
        if GOOGLE_CLIENT_ID:
            # Verify token with audience when client id is provided
            idinfo = id_token.verify_oauth2_token(payload.token, requests.Request(), GOOGLE_CLIENT_ID)
        else:
            # Fallback: verify without audience (less strict)
            idinfo = id_token.verify_oauth2_token(payload.token, requests.Request())
        user_email = idinfo.get("email")
        user_name = idinfo.get("name")
        return JSONResponse(content={"email": user_email, "name": user_name})
    except Exception as e:
        raise HTTPException(status_code=401, detail=str(e))