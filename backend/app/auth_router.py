from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from google.oauth2 import id_token
from google.auth.transport import requests

router = APIRouter()

@router.post("/google-login")
async def google_login(token: str):
    """
    Verify Google login token and return basic user info
    """
    try:
        idinfo = id_token.verify_oauth2_token(token, requests.Request())
        user_email = idinfo.get("email")
        user_name = idinfo.get("name")
        return JSONResponse(content={"email": user_email, "name": user_name})
    except Exception as e:
        raise HTTPException(status_code=401, detail=str(e))